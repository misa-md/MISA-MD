#include <utils/mpi_utils.h>
#include <logs/logs.h>

#include "simulation.h"
#include "hardware_accelerate.hpp"
#include "world_builder.h"
#include "atom_dump.h"

simulation::simulation() : _p_domain(nullptr), _input(nullptr) {
    pConfigVal = &(ConfigParser::getInstance()->configValues);
//    createDomainDecomposition();
//    collision_step = -1;
}

simulation::~simulation() {
//    delete _p_domain; // see finalize method.
    delete _atom;
    delete _integrator;
    delete _pot;

    delete _input; // delete null pointer has no effect.
}

void simulation::createDomainDecomposition() {
    _finalCheckpoint = true;

    //进行区域分解
    kiwi::logs::v(MASTER_PROCESSOR, "domain", "Initializing GlobalDomain decomposition.\n");
    _p_domain = (new Domain(pConfigVal->phaseSpace,
                            pConfigVal->latticeConst,
                            pConfigVal->cutoffRadiusFactor))
            ->decomposition()
            ->createGlobalDomain() // set global box domain.
            ->createSubBoxDomain(); // set local sub-box domain.
    kiwi::logs::v(MASTER_PROCESSOR, "domain", "Initialization done.\n");

//    _numberOfTimesteps = 1;
}

void simulation::createAtoms() {
    _atom = new atom(_p_domain, pConfigVal->latticeConst,
                     pConfigVal->cutoffRadiusFactor, pConfigVal->createSeed);

    if (pConfigVal->createPhaseMode) {  //创建原子坐标、速度信息
        WorldBuilder mWorldBuilder;
        mWorldBuilder.setDomain(_p_domain)
                .setAtomsContainer(_atom)
                .setBoxSize(pConfigVal->phaseSpace[0], pConfigVal->phaseSpace[1], pConfigVal->phaseSpace[2])
                .setRandomSeed(pConfigVal->createSeed)
                .setLatticeConst(pConfigVal->latticeConst)
                .setTset(pConfigVal->createTSet)
                .setAlloyRatio(pConfigVal->alloyRatio)
                .build();
    } else { //读取原子坐标、速度信息
        _input = new input();
        _input->readPhaseSpace(_atom);
    }
    _integrator = new integrator(DEFAULT_TIME_STEP_LENGTH); // time step length.
}

void simulation::prepareForStart() {
    double starttime, stoptime;
    double commtime, computetime, comm;
    _pot = new eam();
    // 读取势函数文件
    if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
        if (!(_pot->parse(pConfigVal->potentialFilename, pConfigVal->potentialFileType))) { // parse potential file.
            abort(1); // todo log, reason
        }
    }

    _pot->eamBCast(kiwi::mpiUtils::own_rank); // BCast Potential
    _pot->interpolateFile(); // interpolation.

    beforeAccelerateRun(_pot); // it runs after atom and boxes creation, but before simulation running.

    starttime = MPI_Wtime();
    _p_domain->exchangeAtomFirst(_atom);
    stoptime = MPI_Wtime();
    commtime = stoptime - starttime;

    _atom->clearForce(); // clear force before running simulation.
    starttime = MPI_Wtime();
    _atom->computeEam(_pot, _p_domain, comm);
    stoptime = MPI_Wtime();
    computetime = stoptime - starttime - comm;
    commtime += comm;

    //_atom->print_force();
    starttime = MPI_Wtime();
    _p_domain->sendForce(_atom);
    stoptime = MPI_Wtime();
    commtime += stoptime - starttime;

    kiwi::logs::i(MASTER_PROCESSOR, "sim", "first step comm time: {}\n", commtime);
    kiwi::logs::i(MASTER_PROCESSOR, "sim", "first step compute time: {}\n", computetime);
}

void simulation::simulate() {
    // start do simulation.
    double starttime, stoptime;
    double commtime = 0, computetime = 0, comm;
    double alltime, allstart, allstop;

    allstart = MPI_Wtime();
    for (_simulation_time_step = 0;
         _simulation_time_step < pConfigVal->timeSteps; _simulation_time_step++) {
        if (_simulation_time_step == pConfigVal->collisionStep) {
            _atom->setv(pConfigVal->collisionLat, pConfigVal->collisionV);
            _p_domain->exchangeInter(_atom);
            _p_domain->borderInter(_atom);
            _p_domain->exchangeAtom(_atom);
            _atom->clearForce();
            _atom->computeEam(_pot, _p_domain, comm);
            _p_domain->sendForce(_atom);
        }
        //先进行求解牛顿运动方程第一步
        _integrator->firststep(_atom);

        //判断是否有粒子跑出晶格点
        kiwi::logs::v(MASTER_PROCESSOR, "simulation", "start deciding atoms:\n");
        _atom->decide();

        //通信ghost区域，交换粒子
        if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
            kiwi::logs::v("simulation", "start ghost communication:\n");
        }
        starttime = MPI_Wtime();
        _p_domain->exchangeInter(_atom);
        _p_domain->borderInter(_atom);
        _p_domain->exchangeAtom(_atom);
        stoptime = MPI_Wtime();
        commtime += stoptime - starttime;

        //计算力
        if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
            kiwi::logs::v("simulation", "start calculating force:\n");
        }
        _atom->clearForce();
        starttime = MPI_Wtime();
        _atom->computeEam(_pot, _p_domain, comm);
        stoptime = MPI_Wtime();
        computetime += stoptime - starttime - comm;
        commtime += comm;

        //发送力
        kiwi::logs::v(MASTER_PROCESSOR, "simulation", "start sending force:\n");
        starttime = MPI_Wtime();
        _p_domain->sendForce(_atom);
        stoptime = MPI_Wtime();
        commtime += stoptime - starttime;
        //求解牛顿运动方程第二步
        _integrator->secondstep(_atom);
    }
    if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
        kiwi::logs::i("simulation", "loop comm time: {}.\n", commtime);
        kiwi::logs::i("simulation", "loop compute time: {}.\n", computetime);
    }
    //输出原子信息
//    if(_simulation_time_step == 10) // todo output atoms every 10 steps.
    output();
    allstop = MPI_Wtime();
    alltime = allstop - allstart;
    if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
        kiwi::logs::i("simulation", "total time:{}.\n", alltime);
    }
}

void simulation::finalize() {
    if (_p_domain != nullptr) {
        delete _p_domain;
        _p_domain = nullptr;
    }
}

void simulation::output() {
    // atom boundary in array.
    _type_lattice_coord begin[DIMENSION] = {
            _p_domain->getSubBoxLatticeCoordLower(0) - _p_domain->getGhostLatticeCoordLower(0),
            _p_domain->getSubBoxLatticeCoordLower(1) - _p_domain->getGhostLatticeCoordLower(1),
            _p_domain->getSubBoxLatticeCoordLower(2) - _p_domain->getGhostLatticeCoordLower(2)};
    _type_lattice_coord end[DIMENSION] = {
            begin[0] + _p_domain->getSubBoxLatticeSize(0),
            begin[1] + _p_domain->getSubBoxLatticeSize(1),
            begin[2] + _p_domain->getSubBoxLatticeSize(2)};
    _type_lattice_size atoms_size = _p_domain->getSubBoxLatticeSize(0) * _p_domain->getSubBoxLatticeSize(1) *
                                    _p_domain->getSubBoxLatticeSize(2);

    AtomDump dump;
    // config dump.
    dump.setDumpFile(pConfigVal->outputDumpFilename)
            .setMode(pConfigVal->outputMode)
            .setBoundary(begin, end, atoms_size);
    dump.dump(_atom);
}

void simulation::abort(int exitcode) {
    MPI_Abort(MPI_COMM_WORLD, exitcode);
}
