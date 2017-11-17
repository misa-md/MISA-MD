#include <iostream>
#include "crystal_md.h"
#include "config.h"
#include "mpi_utils.h"
#include "include/args.hpp"

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
            argvStatus = 1;
        }
        catch (args::ParseError e) {
            cerr << e.what() << endl;
            cerr << parser;
            argvStatus = 2;
        }
        catch (args::ValidationError e) {
            cerr << e.what() << endl;
            cerr << parser;
            argvStatus = 2;
        }

        //check configure
        if (argvStatus == 0) {
            if (conf) {
                config *cp = config::newInstance(args::get(conf));//initial Config
                if (cp->hasError) {
                    cerr << cp->errorMessage;
                    argvStatus = 2;
                }
            } else {
                cerr << "Error: The config file is required" << endl;
                cerr << parser;
                argvStatus = 2;
            }
        }
    } //end if

    // synchronize argvStatus to all processors.
    MPI_Bcast(&argvStatus, 1, MPI_SHORT, MASTER_PROCESSOR, MPI_COMM_WORLD);

    if (argvStatus == 0) { //right argv,and right configure
        cp = config::newInstance();
        MPI_Bcast(cp, sizeof(*cp), MPI_BYTE, MASTER_PROCESSOR, MPI_COMM_WORLD); // synchronize config information
        config::onPostMPICopy(cp);
        //configure check
        return cp->configureCheck();
    } else { //has error or help mode
        return false;
    }
}

bool crystalMD::prepare() {
    simu = new simulation();
    simu->domainDecomposition(); //区域分解
    simu->createBoxedAndAtoms();
    return true; //todo
}

void crystalMD::run() {
    simu->prepareForStart(mpiUtils::ownRank);
    if (mpiUtils::ownRank == MASTER_PROCESSOR)
        std::cout << "Start simulation" << std::endl;
    //开始模拟
    simu->simulate();
}

void crystalMD::destroy() {
    if (mpiUtils::ownRank == MASTER_PROCESSOR)
        std::cout << "finalizing simulation" << std::endl;
    //模拟结束
    simu->finalize();
}

void crystalMD::detach() {
    if (mpiUtils::ownRank == MASTER_PROCESSOR && argvStatus == 0) {
        cout << "app was detached" << endl;
    }
    mpiUtils::finishMPI();
}
