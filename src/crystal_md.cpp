//
// Created by gensh(genshenchu@gmail.com) on 2017/4/19.
//

#include <iostream>
#include <args.hpp>
#include <utils/mpi_utils.h>

#include "crystal_md.h"
#include "arch_env.hpp"

using namespace std;

bool crystalMD::beforeCreate(int argc, char *argv[]) {
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
        return false;
    }
    catch (args::ParseError e) {
        cerr << e.what() << endl;
        cerr << parser;
        return false;
    }
    catch (args::ValidationError e) {
        cerr << e.what() << endl;
        cerr << parser;
        return false;
    }

    // todo add version command.

    if (conf) {
        configFilePath = args::get(conf);
        return true;
    }

    // if no args, print usage.
    std::cerr << parser;
    return false;
}

void crystalMD::onCreate() {
    Config *pConfig;
    if (kiwi::mpiUtils::ownRank == MASTER_PROCESSOR) {
        std::cout << "mpi env was initialed." << std::endl;
        // initial config Obj, then read and resolve config file.
        pConfig = Config::newInstance(configFilePath); // todo config file from argv.
        if (pConfig->hasError) {
            std::cerr << "[Error] " << pConfig->errorMessage << std::endl;
            this->abort(2);
        }
    } else {
        // just initial a empty config Obj.
        pConfig = Config::getInstance();
    }
    pConfig->sync(); // sync config data to other processors from master processor.
#ifdef DEV_MODE
// print configure.
    if (kiwi::mpiUtils::ownRank != MASTER_PROCESSOR) {
        cout << pConfig->configValues;
    }
#endif
    archEnvInit(); // initialize architectures environment.
}

bool crystalMD::prepare() {
    pSimulation = new simulation();
    pSimulation->domainDecomposition(); //区域分解
    pSimulation->createBoxedAndAtoms();
    return true;
}

void crystalMD::onStart() {
    pSimulation->prepareForStart(kiwi::mpiUtils::ownRank);
    if (kiwi::mpiUtils::ownRank == MASTER_PROCESSOR)
        std::cout << "Start simulation" << std::endl;
    //开始模拟
    pSimulation->simulate();
}

void crystalMD::onFinish() {
    if (kiwi::mpiUtils::ownRank == MASTER_PROCESSOR)
        std::cout << "finalizing simulation" << std::endl;
    //模拟结束
    pSimulation->finalize();
}

void crystalMD::beforeDestroy() {
    if (kiwi::mpiUtils::ownRank == MASTER_PROCESSOR) {
        cout << "app was detached" << endl;
    }
    archEnvFinalize(); // clean architectures environment.
}

void crystalMD::onDestroy() {}
