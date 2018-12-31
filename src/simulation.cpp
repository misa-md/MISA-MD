#include <utils/mpi_utils.h>
#include <logs/logs.h>
#include <eam.h>
#include <parser/setfl_parser.h>

#include "simulation.h"
#include "utils/mpi_domain.h"
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
    delete _atom;
    delete _newton_motion;
    delete _pot;

    delete _input; // delete null pointer has no effect.
}

void simulation::createDomainDecomposition() {
    _finalCheckpoint = true;

    //进行区域分解
    kiwi::logs::v(MASTER_PROCESSOR, "domain", "Initializing GlobalDomain decomposition.\n");
    MPI_Comm new_comm;
    _p_domain = Domain::Builder()
            .setComm(MPIDomain::sim_processor, &new_comm)
            .setPhaseSpace(pConfigVal->phaseSpace)
            .setCutoffRadius(pConfigVal->cutoffRadiusFactor)
            .setLatticeConst(pConfigVal->latticeConst)
            .build();
    kiwi::mpiUtils::onGlobalCommChanged(new_comm); // set new domain.
    MPIDomain::sim_processor = kiwi::mpiUtils::global_process;

    kiwi::logs::v(MASTER_PROCESSOR, "domain", "Initialization done.\n");
}

void simulation::createAtoms() {
    _atom = new atom(_p_domain, pConfigVal->latticeConst, pConfigVal->cutoffRadiusFactor);
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

    // todo file type funl support. pConfigVal->potentialFileType
    // 读取势函数文件
    //atom_type::_type_atom_types eles = 0;
    if (MPIDomain::sim_processor.own_rank == MASTER_PROCESSOR) {
        char tmp[4096];
        sprintf(tmp, "%s", pConfigVal->potentialFilename.c_str());

        FILE *pot_file = fopen(tmp, "r");
        if (pot_file == nullptr) { // todo open too many in md.
            kiwi::logs::e("pot", "open potential file {} failed.\n", pConfigVal->potentialFilename);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return;
        }
        // new parser
        SetflParser *parser = new SetflParser(pConfigVal->potentialFilename); // todo delete (vector)
        parser->parseHeader(); // elements count got. // todo parsing error.
        // eles = parser->getEles(); // elements values on non-root processors are 0.
        _pot = eam::newInstance(parser->getEles(),
                                MASTER_PROCESSOR,
                                MPIDomain::sim_processor.own_rank,
                                MPIDomain::sim_processor.comm);
        // read data
        parser->parseBody(_pot); // todo parsing error.
        parser->done();
/*      BCast type list, the lists size is the same as eam potential elements size.
        parser->type_lists.sync(MASTER_PROCESSOR,
                                MPIDomain::sim_processor.own_rank,
                                MPIDomain::sim_processor.comm,
                                _pot->geEles());
*/
    } else {
        _pot = eam::newInstance(0,
                                MASTER_PROCESSOR,
                                MPIDomain::sim_processor.own_rank,
                                MPIDomain::sim_processor.comm);
    }
    // BCast Potential
    _pot->eamBCast(MASTER_PROCESSOR,
                   MPIDomain::sim_processor.own_rank,
                   MPIDomain::sim_processor.comm);
    _pot->interpolateFile(); // interpolation.

    beforeAccelerateRun(_pot); // it runs after atom and boxes creation, but before simulation running.

    starttime = MPI_Wtime();
    // todo make _cut_lattice a member of class AtomList
    _atom->getAtomList()->exchangeAtomFirst(_p_domain, _atom->getCutLattice());
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
    _atom->sendForce();
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
            if (!pConfigVal->originDumpPath.empty()) {
                output(_simulation_time_step, true); // dump atoms
            }
            _atom->setv(pConfigVal->collisionLat, pConfigVal->direction, pConfigVal->pkaEnergy);
            _atom->getInterList()->exchangeInter(_p_domain);
            _atom->getInterList()->borderInter(_p_domain);
            _atom->getAtomList()->exchangeAtom(_p_domain);
            _atom->clearForce();
            _atom->computeEam(_pot, _p_domain, comm);
            _atom->sendForce();
        }
        //先进行求解牛顿运动方程第一步
        _newton_motion->firststep(_atom->getAtomList(), _atom->getInterList());

        //判断是否有粒子跑出晶格点
        _atom->decide();

        //通信ghost区域，交换粒子
        starttime = MPI_Wtime();
        _atom->getInterList()->exchangeInter(_p_domain);
        _atom->getInterList()->borderInter(_p_domain);
        _atom->getAtomList()->exchangeAtom(_p_domain);
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
        _atom->sendForce();
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

void simulation::output(size_t time_step, bool before_collision) {
    // atom boundary in array.
    _type_lattice_coord begin[DIMENSION] = {
            _p_domain->lattice_coord_sub_box_region.x_low - _p_domain->lattice_coord_ghost_region.x_low,
            _p_domain->lattice_coord_sub_box_region.y_low - _p_domain->lattice_coord_ghost_region.y_low,
            _p_domain->lattice_coord_sub_box_region.z_low - _p_domain->lattice_coord_ghost_region.z_low};
    _type_lattice_coord end[DIMENSION] = {
            begin[0] + _p_domain->lattice_size_sub_box[0],
            begin[1] + _p_domain->lattice_size_sub_box[1],
            begin[2] + _p_domain->lattice_size_sub_box[2]};
    _type_lattice_size atoms_size = _p_domain->lattice_size_sub_box[0] * _p_domain->lattice_size_sub_box[1] *
                                    _p_domain->lattice_size_sub_box[2];
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
        std::string filename;
        if (before_collision) {
            filename = pConfigVal->originDumpPath; // todo pass file name from func output parameters.
        } else {
            filename = fmt::format(pConfigVal->atomsDumpFilePath, time_step);
        }
        // pointer to the atom dump class for outputting atoms information.
        auto *dumpInstance = new AtomDump(pConfigVal->atomsDumpMode, filename,
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
