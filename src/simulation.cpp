#include <iostream>
#include "simulation.h"
#include "domain.h"
#include "domaindecomposition.h"
#include "atom.h"
#include "mpi_utils.h"
#include "config.h"

simulation::simulation() : _domaindecomposition(nullptr), _input(nullptr) {
//    domainDecomposition();

    //collision_step = -1;
    cp = config::newInstance();
}

simulation::~simulation() {
    delete _domain;
    delete _atom;
    delete _integrator;
    delete _pot;
    if (_input != nullptr)
        delete _input;
    if (_createatom != nullptr)
        delete _createatom;
}

void simulation::domainDecomposition() {
//    _domaindecomposition = NULL;
    _domain = nullptr;
    _finalCheckpoint = true;

    //进行区域分解
    if (mpiUtils::ownRank == MASTER_PROCESSOR)
        std::cout << "Initializing domain decomposition ... " << std::endl;
    _domaindecomposition = new domaindecomposition();
    if (mpiUtils::ownRank == MASTER_PROCESSOR)
        std::cout << "Initialization done" << std::endl;

    if (mpiUtils::ownRank == MASTER_PROCESSOR)
        std::cout << "Constructing domain ..." << std::endl;
    _domain = new domain(mpiUtils::ownRank);
    if (mpiUtils::ownRank == MASTER_PROCESSOR)
        std::cout << "Domain construction done." << std::endl;

//    /*
//     * 初始化参数
//     */
//    _numberOfTimesteps = 1;
//
}

void simulation::createBoxedAndAtoms() {
    double mass, ghostLength;

    double boxlo[3], boxhi[3], globalLength[3];
    boxlo[0] = boxlo[1] = boxlo[2] = 0;
    globalLength[0] = boxhi[0] = cp->phaseSpace[0] * cp->latticeConst; //box_x个单位长度(单位长度即latticeconst)
    globalLength[1] = boxhi[1] = cp->phaseSpace[1] * cp->latticeConst;
    globalLength[2] = boxhi[2] = cp->phaseSpace[2] * cp->latticeConst;

    _domain->setGlobalLength(0, globalLength[0]);
    _domain->setGlobalLength(1, globalLength[1]);
    _domain->setGlobalLength(2, globalLength[2]);

    double bBoxMin[3], bBoxMax[3];
    for (int i = 0; i < 3; i++) {
        bBoxMin[i] = _domaindecomposition->getBoundingBoxMin(i, _domain);
        bBoxMax[i] = _domaindecomposition->getBoundingBoxMax(i, _domain);
    }
    ghostLength = cp->cutoffRadius;
    _atom = new atom(boxlo, boxhi, globalLength, bBoxMin, bBoxMax, ghostLength,
                     cp->latticeConst, cp->cutoffRadius, cp->createSeed);
    mass = 55.845;

    if (cp->createPhaseMode) {  //创建原子坐标、速度信息
        _createatom = new createatom(cp->createTSet);
        _createatom->createphasespace(_atom, mass, cp->phaseSpace[0], cp->phaseSpace[1], cp->phaseSpace[2]);
    } else { //读取原子坐标、速度信息
        _input = new input();
        _input->readPhaseSpace(_atom);
    }
    _integrator = new integrator(0.001);
}

void simulation::prepareForStart(int rank) {
    double starttime, stoptime;
    double commtime, computetime, comm;
    _pot = new eam();
    //读取势函数文件
    if (rank == 0) {
        initEamPotential(cp->potentialFileType);
    }
    eamBCastPotential(rank);
    eamPotentialInterpolate();

    starttime = MPI_Wtime();
    _domaindecomposition->exchangeAtomfirst(_atom, _domain);
    stoptime = MPI_Wtime();
    commtime = stoptime - starttime;
    starttime = MPI_Wtime();
    _atom->computeEam(_pot, _domaindecomposition, comm);
    stoptime = MPI_Wtime();
    computetime = stoptime - starttime - comm;
    commtime += comm;
    //_atom->print_force();
    starttime = MPI_Wtime();
    _domaindecomposition->sendforce(_atom);
    stoptime = MPI_Wtime();
    commtime += stoptime - starttime;
    if (mpiUtils::ownRank == MASTER_PROCESSOR) {
        printf("first step comm time: %lf\n", commtime);
        printf("first step compute time: %lf\n", computetime);
    }
}

void simulation::simulate() {
    //开始进行模拟
    double starttime, stoptime;
    double commtime = 0, computetime = 0, comm;
    int nflag;
    for (_simstep = 0; _simstep < cp->timeSteps; _simstep++) {
        if (_simstep == cp->collisionSteps) {
            _atom->setv(cp->collisionLat, cp->collisionV);
            _domaindecomposition->exchangeInter(_atom, _domain);
            _domaindecomposition->borderInter(_atom, _domain);
            _domaindecomposition->exchangeAtom(_atom, _domain);
            _atom->clear_force();
            _atom->computeEam(_pot, _domaindecomposition, comm);
            _domaindecomposition->sendforce(_atom);
        }
        //先进行求解牛顿运动方程第一步
        _integrator->firststep(_atom);
        //判断是否有粒子跑出晶格点
        nflag = _atom->decide();

        //通信ghost区域，交换粒子
        starttime = MPI_Wtime();
        _domaindecomposition->exchangeInter(_atom, _domain);
        _domaindecomposition->borderInter(_atom, _domain);
        _domaindecomposition->exchangeAtom(_atom, _domain);
        stoptime = MPI_Wtime();
        commtime += stoptime - starttime;

        //计算力
        _atom->clear_force();
        starttime = MPI_Wtime();
        _atom->computeEam(_pot, _domaindecomposition, comm);
        stoptime = MPI_Wtime();
        computetime += stoptime - starttime - comm;
        commtime += comm;
        //发送力
        starttime = MPI_Wtime();
        _domaindecomposition->sendforce(_atom);
        stoptime = MPI_Wtime();
        commtime += stoptime - starttime;
        //求解牛顿运动方程第二步
        _integrator->secondstep(_atom);
    }
    if (mpiUtils::ownRank == MASTER_PROCESSOR) {
        printf("loop comm time: %lf\n", commtime);
        printf("loop compute time: %lf\n", computetime);
    }
    //输出原子信息
    output();
}

void simulation::finalize() {
    if (_domaindecomposition != nullptr) {
        delete _domaindecomposition;
        _domaindecomposition = nullptr;
    }
}

//读取势函数文件，共分为两种势函数文件格式：
//"funcfl"和"setfl"
void simulation::initEamPotential(string file_type) {
    if (file_type == "funcfl") {
        char tmp[4096];
        sprintf(tmp, "%s", cp->potentialFilename.c_str());
        FILE *potFile = fopen(tmp, "r");
        if (potFile == nullptr) {
            std::cerr << "file not found" << std::endl;
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
        sprintf(tmp, "%s", cp->potentialFilename.c_str());

        FILE *potFile = fopen(tmp, "r");
        if (potFile == NULL) {
            std::cerr << "file not found!" << std::endl;
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
        if (strtok(copy, " \t\n\r\f") == NULL) {
            n = 0;
        } else {
            n = 1;
            while (strtok(NULL, " \t\n\r\f")) n++;
        }
        int nwords = n;
        delete[] copy;
        if (nwords != nElems + 1)
            std::cerr << "Incorrect element names in EAM potential file!" << endl;

        char **words = new char *[nElems + 1];
        nwords = 0;
        strtok(tmp, " \t\n\r\f");
        while ((words[nwords++] = strtok(NULL, " \t\n\r\f"))) continue;

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
    int ownrank = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &ownrank);
    _atom->print_atom(ownrank);
}

void simulation::exit(int exitcode) {
    MPI_Abort(MPI_COMM_WORLD, exitcode);
}
