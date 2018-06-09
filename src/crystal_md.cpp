//
// Created by gensh(genshenchu@gmail.com) on 2017/4/19.
//

#include <args.hpp>
#include <utils/mpi_utils.h>
#include <logs/logs.h>

#include "crystal_md.h"
#include "arch_env.hpp"

bool crystalMD::beforeCreate(int argc, char *argv[]) {
    // parser arguments
    // see https://github.com/Taywee/args for using args.
    args::ArgumentParser parser("This is CrystalMD program.", "authors:BaiHe.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> conf(parser, "conf", "The configure file", {'c', "conf"});
    try {
        parser.ParseCLI(argc, (const char *const *) argv);
    }
    catch (args::Help) {
        std::cout << parser;
        return false;
    }
    catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return false;
    }
    catch (args::ValidationError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
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
    ConfigParser *pConfig;
    if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
        kiwi::logs::s("env", "mpi env is initialized.\n");
        // initial config Obj, then read and resolve config file.
        pConfig = ConfigParser::newInstance(configFilePath); // todo config file from argv.
        if (pConfig->hasError) {
            kiwi::logs::e("config", "{0}\n", pConfig->errorMessage);
            this->abort(2);
        }
    } else {
        // just initial a empty config Obj.
        pConfig = ConfigParser::getInstance();
    }
    pConfig->sync(); // sync config data to other processors from master processor.
#ifdef DEV_MODE
// print configure.
    if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
        std::cout << pConfig->configValues;
    }
#endif

    // prepare logs.
    if (pConfig->configValues.logs_mode == LOGS_MODE_CONSOLE) {
        kiwi::logs::setCorlorFul(true);
    } else if (pConfig->configValues.logs_mode == LOGS_MODE_FILE) {
        kiwi::logs::setLogFile(pConfig->configValues.logs_filename);
    }

    archEnvInit(); // initialize architectures environment.
}

bool crystalMD::prepare() {
    pSimulation = new simulation();
    pSimulation->createDomainDecomposition(); // 区域分解
    pSimulation->createAtoms();
    return true;
}

void crystalMD::onStart() {
    pSimulation->prepareForStart();
    kiwi::logs::v(MASTER_PROCESSOR, "simulation", "Start simulation.\n");
    pSimulation->simulate(); // start simulation.
}

void crystalMD::onFinish() {
    kiwi::logs::s(MASTER_PROCESSOR, "simulation", "finalizing simulation\n");
    //模拟结束
    pSimulation->finalize();
}

void crystalMD::beforeDestroy() {
    kiwi::logs::v(MASTER_PROCESSOR, "app", "app was detached.\n");
    archEnvFinalize(); // clean architectures environment.
}

void crystalMD::onDestroy() {}
