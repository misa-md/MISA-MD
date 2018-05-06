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
//    delete p_domain;
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
    const double mass = 55.845;

    if (pConfigVal->createPhaseMode) {  //创建原子坐标、速度信息
        WorldBuilder mWorldBuilder;
        mWorldBuilder.setDomain(_p_domain)
                .setAtomsContainer(_atom)
                .setBoxSize(pConfigVal->phaseSpace[0], pConfigVal->phaseSpace[1], pConfigVal->phaseSpace[2])
                .setRandomSeed(pConfigVal->createSeed)
                .setLatticeConst(pConfigVal->latticeConst)
                .setTset(pConfigVal->createTSet)
                .setMass(mass)
                .build();
    } else { //读取原子坐标、速度信息
        _input = new input();
        _input->readPhaseSpace(_atom);
    }
    _integrator = new integrator(0.001); // time step width.
}

void simulation::prepareForStart(int rank) {
    double starttime, stoptime;
    double commtime, computetime, comm;
    _pot = new eam();
    //读取势函数文件
    if (rank == 0) {
        initEamPotential(pConfigVal->potentialFileType);
    }
    eamBCastPotential(rank);
    eamPotentialInterpolate();

    beforeAccelerateRun(_pot); // it runs after atom and boxes creation, but before simulation running.

    starttime = MPI_Wtime();
    _p_domain->exchangeAtomfirst(_atom);
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
    if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
        kiwi::logs::i("sim", "first step comm time: {}\n", commtime);
        kiwi::logs::i("sim", "first step compute time: {}\n", computetime);
    }
}

void simulation::simulate() {
    //开始进行模拟
    double starttime, stoptime;
    double commtime = 0, computetime = 0, comm;
    double alltime, allstart, allstop;
    int nflag;

    allstart = MPI_Wtime();
    for (_simulation_time_step = 0;
         _simulation_time_step < pConfigVal->timeSteps; _simulation_time_step++) {
        if (_simulation_time_step == pConfigVal->collisionSteps) {
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
        if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
            kiwi::logs::v("simulation", "start deciding atoms:\n");
        }
        nflag = _atom->decide();

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
        if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
            kiwi::logs::v("simulation", "start sending force:\n");
        }
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

//读取势函数文件，共分为两种势函数文件格式：
//"funcfl"和"setfl"
void simulation::initEamPotential(string file_type) {
    if (file_type == "funcfl") {
        char tmp[4096];
        sprintf(tmp, "%s", pConfigVal->potentialFilename.c_str());
        FILE *potFile = fopen(tmp, "r");
        if (potFile == nullptr) {
            kiwi::logs::e("eam", "file {} not found.\n", pConfigVal->potentialFilename);
            exit(1);
        }

        // 第一行
        fgets(tmp, sizeof(tmp), potFile);
        char name[3];
        sscanf(tmp, "%s", name);
        _pot->setname(name);

        // 第二行
        int nAtomic;
        double mass, lat;
        char latticeType[8];
        fgets(tmp, sizeof(tmp), potFile);
        sscanf(tmp, "%d %le %le %s", &nAtomic, &mass, &lat, latticeType);
        _pot->setatomicNo(nAtomic); // 原子序号
        _pot->setlat(lat); // 晶格常数
        _pot->setmass(0, mass); // 质量.
        _pot->setlatticeType(latticeType); //晶格类型

        // 第三行
        int nRho, nR;
        double dRho, dR, cutoff;
        fgets(tmp, sizeof(tmp), potFile);
        sscanf(tmp, "%d %le %d %le %le", &nRho, &dRho, &nR, &dR, &cutoff);
        _pot->setcutoff(cutoff); //截断半径
        double x0 = 0.0;

        // 申请读取数据的空间
        int bufSize = max(nRho, nR);
        double *buf = new double[bufSize];

        // 读取嵌入能表
        for (int ii = 0; ii < nRho; ++ii)
            fscanf(potFile, "%lg", buf + ii);
        _pot->initf(0, nRho, x0, dRho, buf); //通过读取势文件的数据建立table

        // 读取对势表
        for (int ii = 0; ii < nR; ++ii)
            fscanf(potFile, "%lg", buf + ii);
        double r;
        for (int ii = 1; ii < nR; ++ii) {
            r = x0 + ii * dR;
            buf[ii] *= buf[ii] / r;
            buf[ii] *= hartreeToEv * bohrToAngs;
        }
        buf[0] = buf[1] + (buf[1] - buf[2]);
        _pot->initphi(0, nR, x0, dR, buf);

        // 读取电子云密度表
        for (int ii = 0; ii < nR; ++ii)
            fscanf(potFile, "%lg", buf + ii);
        _pot->initrho(0, nR, x0, dR, buf);

        delete[] buf;
    } else if (string(file_type) == string("setfl")) {
        char tmp[4096];
        sprintf(tmp, "%s", pConfigVal->potentialFilename.c_str());

        FILE *potFile = fopen(tmp, "r");
        if (potFile == nullptr) {
            kiwi::logs::e("eam", "file {} not found.\n", pConfigVal->potentialFilename);
            exit(1);
        }

        // 前三行为注释
        fgets(tmp, sizeof(tmp), potFile);
        fgets(tmp, sizeof(tmp), potFile);
        fgets(tmp, sizeof(tmp), potFile);

        // 第四行
        fgets(tmp, sizeof(tmp), potFile);
        int nElems;
        sscanf(tmp, "%d", &nElems); //原子类型个数

        _pot->init(nElems);//对势函数进行初始化，从文件中读入原子类型个数后。

        char *copy;
        copy = new char[strlen(tmp) + 1];
        strcpy(copy, tmp);
        char *ptr;
        if ((ptr = strchr(copy, '#'))) *ptr = '\0';
        int n;
        if (strtok(copy, " \t\n\r\f") == nullptr) {
            n = 0;
        } else {
            n = 1;
            while (strtok(nullptr, " \t\n\r\f")) n++;
        }
        int nwords = n;
        delete[] copy;
        if (nwords != nElems + 1) {
            kiwi::logs::e("eam", "Incorrect element names in EAM potential file!");
            // todo MPI abort.
        }

        char **words = new char *[nElems + 1];
        nwords = 0;
        strtok(tmp, " \t\n\r\f");
        while ((words[nwords++] = strtok(nullptr, " \t\n\r\f"))) continue;

        delete[] words;
        // 第五行
        int nRho, nR;
        double dRho, dR, cutoff;
        // 所有原子使用同一个截断半径
        fgets(tmp, sizeof(tmp), potFile);
        sscanf(tmp, "%d %le %d %le %le", &nRho, &dRho, &nR, &dR, &cutoff);
        _pot->setcutoff(cutoff);

        // 申请读取数据空间
        int bufSize = max(nRho, nR);
        double *buf = new double[bufSize];
        double x0 = 0.0;
        // 每种原子信息
        for (int i = 0; i < nElems; i++) {
            fgets(tmp, sizeof(tmp), potFile);
            int nAtomic;
            double mass, lat;
            char latticeType[8];
            sscanf(tmp, "%d %le %le %s", &nAtomic, &mass, &lat, latticeType);

            _pot->setmass(i, mass);  // 原子质量

            // 读取嵌入能表
            grab(potFile, nRho, buf);
            _pot->initf(i, nRho, x0, dRho, buf);

            // 读取电子云密度表
            grab(potFile, nR, buf);
            _pot->initrho(i, nR, x0, dR, buf);
        }

        //读取对势表
        int i, j, k = 0;
        for (i = 0; i < nElems; i++) {
            for (j = 0; j <= i; j++) {
                grab(potFile, nR, buf);
                _pot->initphi(k++, nR, x0, dR, buf);
            }
        }
        delete[] buf;
    }
}

void simulation::grab(FILE *fptr, int n, double *list) {
    char *ptr;
    char line[1024];

    int i = 0;
    while (i < n) {
        fgets(line, 1024, fptr);
        ptr = strtok(line, " \t\n\r\f");
        list[i++] = atof(ptr);
        while ((ptr = strtok(nullptr, " \t\n\r\f"))) list[i++] = atof(ptr);
    }
}

void simulation::eamBCastPotential(int rank) {
    _pot->eamBcast(rank);
}

void simulation::eamPotentialInterpolate() {
    _pot->interpolatefile();
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

void simulation::exit(int exitcode) {
    MPI_Abort(MPI_COMM_WORLD, exitcode);
}
