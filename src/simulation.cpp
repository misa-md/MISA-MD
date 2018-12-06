#include <utils/mpi_utils.h>
#include <logs/logs.h>

#include "simulation.h"
#include "utils/mpi_domain.h"
#include "potential/eam_parser.h"
#include "hardware_accelerate.hpp"
#include "world_builder.h"
#include "atom_dump.h"

simulation::simulation() : _p_domain(nullptr), _atom(nullptr),
                           _newton_motion(nullptr), _input(nullptr), _pot(nullptr) {
    pConfigVal = &(ConfigParser::getInstance()->configValues);
//    createDomainDecomposition();
//    collision_step = -1;
}

simulation::~simulation() {
//    delete _p_domain; // see finalize method.
    delete _atom;
    delete _newton_motion;
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
    _atom->calculateNeighbourIndices(); // establish index offset for neighbour.

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
    _newton_motion = new NewtonMotion(pConfigVal->timeStepLength); // time step length.
}

void simulation::prepareForStart() {
    double starttime, stoptime;
    double commtime, computetime, comm;
    _pot = new eam();
    EamParser parser(pConfigVal->potentialFilename, pConfigVal->potentialFileType);
    // 读取势函数文件
    if (MPIDomain::sim_processor.own_rank == MASTER_PROCESSOR) {
        if (!(parser.parse(_pot))) { // parse potential file.
            abort(1); // todo log, reason
        }
    }

    _pot->eamBCast(MPIDomain::sim_processor.own_rank); // BCast Potential
    _pot->interpolateFile(); // interpolation.

    beforeAccelerateRun(_pot); // it runs after atom and boxes creation, but before simulation running.

    starttime = MPI_Wtime();
    _p_domain->exchangeAtomFirst(_atom);
    // fixme those code does not fit [read atom mode], because in [create atom mode], inter is empty at first step;
    // so borderInter is not get called.
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
    for (_simulation_time_step = 0; _simulation_time_step < pConfigVal->timeSteps; _simulation_time_step++) {
        kiwi::logs::s(MASTER_PROCESSOR, "simulation", "simulating steps: {}/{}\r",
                      _simulation_time_step + 1, pConfigVal->timeSteps);

        if (_simulation_time_step == pConfigVal->collisionStep) {
            _atom->setv(pConfigVal->collisionLat, pConfigVal->direction, pConfigVal->pkaEnergy);
            _p_domain->exchangeInter(_atom);
            _p_domain->borderInter(_atom);
            _p_domain->exchangeAtom(_atom);
            _atom->clearForce();
            _atom->computeEam(_pot, _p_domain, comm);
            _p_domain->sendForce(_atom);
        }
        //先进行求解牛顿运动方程第一步
        _newton_motion->firststep(_atom->getAtomList(), _atom->getInterList());

        //判断是否有粒子跑出晶格点
        _atom->decide();

        //通信ghost区域，交换粒子
        starttime = MPI_Wtime();
        _p_domain->exchangeInter(_atom);
        _p_domain->borderInter(_atom);
        _p_domain->exchangeAtom(_atom);
        stoptime = MPI_Wtime();
        commtime += stoptime - starttime;

        //计算力
        _atom->clearForce();
        starttime = MPI_Wtime();
        _atom->computeEam(_pot, _p_domain, comm);
        stoptime = MPI_Wtime();
        computetime += stoptime - starttime - comm;
        commtime += comm;

        //发送力
        starttime = MPI_Wtime();
        _p_domain->sendForce(_atom);
        stoptime = MPI_Wtime();
        commtime += stoptime - starttime;
        //求解牛顿运动方程第二步
        _newton_motion->secondstep(_atom->getAtomList(), _atom->getInterList());

        //输出原子信息
        if ((_simulation_time_step + 1) % pConfigVal->atomsDumpInterval == 0) {// todo output atoms every 10 steps.
            output(_simulation_time_step + 1);
        }
    }
    if (MPIDomain::sim_processor.own_rank == MASTER_PROCESSOR) {
        kiwi::logs::i("simulation", "loop comm time: {}.\n", commtime);
        kiwi::logs::i("simulation", "loop compute time: {}.\n", computetime);
    }
    allstop = MPI_Wtime();
    alltime = allstop - allstart;
    if (MPIDomain::sim_processor.own_rank == MASTER_PROCESSOR) {
        kiwi::logs::i("simulation", "total time:{}.\n", alltime);
    }
}

void simulation::finalize() {
    if (_p_domain != nullptr) {
        delete _p_domain;
        _p_domain = nullptr;
    }
}

void simulation::output(size_t time_step) {
    // atom boundary in array.
    _type_lattice_coord begin[DIMENSION] = {
            _p_domain->getGlobalSubBoxLatticeCoordLower(0) - _p_domain->getGlobalGhostLatticeCoordLower(0),
            _p_domain->getGlobalSubBoxLatticeCoordLower(1) - _p_domain->getGlobalGhostLatticeCoordLower(1),
            _p_domain->getGlobalSubBoxLatticeCoordLower(2) - _p_domain->getGlobalGhostLatticeCoordLower(2)};
    _type_lattice_coord end[DIMENSION] = {
            begin[0] + _p_domain->getSubBoxLatticeSize(0),
            begin[1] + _p_domain->getSubBoxLatticeSize(1),
            begin[2] + _p_domain->getSubBoxLatticeSize(2)};
    _type_lattice_size atoms_size = _p_domain->getSubBoxLatticeSize(0) * _p_domain->getSubBoxLatticeSize(1) *
                                    _p_domain->getSubBoxLatticeSize(2);
    double start = 0, stop = 0;
    static double totalDumpTime = 0;

    start = MPI_Wtime();
    if (!pConfigVal->outByFrame) {
        static AtomDump *dumpInstance = nullptr; // pointer used for non-by-frame dumping.
        if (dumpInstance == nullptr) { // initialize atomDump if it is not initialized.
            dumpInstance = new AtomDump(pConfigVal->atomsDumpMode, pConfigVal->atomsDumpFilePath,
                                        begin, end, atoms_size); // atoms dump.
            // fixme Attempting to use an MPI routine after finalizing MPICH.
        }
        dumpInstance->dump(_atom->getAtomList(), _atom->getInterList(), time_step);
        if (time_step + pConfigVal->atomsDumpInterval > pConfigVal->timeSteps) { // the last time of dumping.
            dumpInstance->writeDumpHeader();
            delete dumpInstance;
        }
    } else {
        std::string filename = fmt::format(pConfigVal->atomsDumpFilePath, time_step);
        // pointer to the atom dump class for outputting atoms information.
        AtomDump *dumpInstance = new AtomDump(pConfigVal->atomsDumpMode, filename,
                                              begin, end, atoms_size);
        dumpInstance->dump(_atom->getAtomList(), _atom->getInterList(), time_step);
        dumpInstance->writeDumpHeader();
        delete dumpInstance;
    }
    stop = MPI_Wtime();
    totalDumpTime += (stop - start);

    // log dumping time.
    if (time_step + pConfigVal->atomsDumpInterval > pConfigVal->timeSteps &&
        MPIDomain::sim_processor.own_rank == MASTER_PROCESSOR) {
        if (pConfigVal->atomsDumpMode == OUTPUT_COPY_MODE) {
            kiwi::logs::i("dump", "time of dumping atoms in copy mode:{}.\n", totalDumpTime);
        } else if (pConfigVal->atomsDumpMode == OUTPUT_DIRECT_MODE) {
            kiwi::logs::i("dump", "time of dumping atoms in direct mode:{}.\n", totalDumpTime);
        }
    }
}

void simulation::abort(int exitcode) {
    MPI_Abort(MPI_COMM_WORLD, exitcode);
}
