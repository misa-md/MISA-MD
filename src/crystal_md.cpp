#include <iostream>
#include "crystal_md.h"
#include "config.h"
#include "mpi_utils.h"
#include "include/args.hpp"
#include "arch_env.hpp"

//
// Created by gensh(genshenchu@gmail.com) on 2017/4/19.
//

using namespace std;

crystalMD::crystalMD(int argc, char **argv) : argc(argc), argv(argv) {}

bool crystalMD::initialize() {
    mpiUtils::initialMPI();

    if (mpiUtils::ownRank == MASTER_PROCESSOR) {
        // parser arguments
        // see https://github.com/Taywee/args for using args.
        args::ArgumentParser parser("This is CrystalMD program.", "authors:BaiHe.");
        args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
        args::ValueFlag<string> conf(parser, "conf", "The configure file", {'c', "conf"});
        try {
            parser.ParseCLI(argc, (const char *const *) argv);
        }
        catch (args::Help) {
            cout << parser;
            mArgvStatus = 1;
        }
        catch (args::ParseError e) {
            cerr << e.what() << endl;
            cerr << parser;
            mArgvStatus = 2;
        }
        catch (args::ValidationError e) {
            cerr << e.what() << endl;
            cerr << parser;
            mArgvStatus = 2;
        }

        //check configure
        if (mArgvStatus == 0) {
            if (conf) {
                config *cp = config::newInstance(args::get(conf));//initial Config
                if (cp->hasError) {
                    cerr << cp->errorMessage;
                    mArgvStatus = 2;
                }
            } else {
                cerr << "Error: The config file is required" << endl;
                cerr << parser;
                mArgvStatus = 2;
            }
        }
    } //end if

    // synchronize mArgvStatus to all processors.
    MPI_Bcast(&mArgvStatus, 1, MPI_SHORT, MASTER_PROCESSOR, MPI_COMM_WORLD);

    if (mArgvStatus == 0) { //right argv,and right configure
        pConfig = config::newInstance();
        MPI_Bcast(pConfig, sizeof(*pConfig), MPI_BYTE, MASTER_PROCESSOR, MPI_COMM_WORLD); // synchronize config information
        config::onPostMPICopy(pConfig);
        //configure check
        if (pConfig->configureCheck()) { // configure check passed.
            return this->runtimeEnvInitialize();
        } else {
            return false;
        }
    } else { //has error or help mode
        return false;
    }
}

/**
 * initial runtime
 */
bool crystalMD::runtimeEnvInitialize() {
    archEnvInit(); // architectures environment initialize.
    return true;
}

bool crystalMD::prepare() {
    pSimulation = new simulation();
    pSimulation->domainDecomposition(); //区域分解
    pSimulation->createBoxedAndAtoms();
    return true; //todo
}

void crystalMD::run() {
    pSimulation->prepareForStart(mpiUtils::ownRank);
    if (mpiUtils::ownRank == MASTER_PROCESSOR)
        std::cout << "Start simulation" << std::endl;
    //开始模拟
    pSimulation->simulate();
}

void crystalMD::destroy() {
    if (mpiUtils::ownRank == MASTER_PROCESSOR)
        std::cout << "finalizing simulation" << std::endl;
    //模拟结束
    pSimulation->finalize();
}

void crystalMD::detach() {
    if (mpiUtils::ownRank == MASTER_PROCESSOR && mArgvStatus == 0) {
        cout << "app was detached" << endl;
    }
    archEnvFinalize(); // clean architectures environment.
    mpiUtils::finishMPI();
}
