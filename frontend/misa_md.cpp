//
// Created by gensh(genshenchu@gmail.com) on 2017/4/19.
//

#include <args.hpp>
#include <utils/mpi_utils.h>
#include <logs/logs.h>
#include <utils/mpi_data_types.h>

#include "misa_md.h"
#include "utils/mpi_domain.h"
#include "arch/arch_env.hpp"
#include "device.h"
#include "frontend_config.h"
#include "md_simulation.h"

bool MISAMD::beforeCreate(int argc, char *argv[]) {
    // parser arguments
    // see https://github.com/Taywee/args for using args.
    args::ArgumentParser parser("This is MISA-MD program.", "authors:BaiHe.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> conf(parser, "conf", "The configure file", {'c', "conf"});
    args::Flag version(parser, "version", "show version number", {'v', "version"});
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

    if (version) {
        std::cout << "MISA-MD version " << MD_VERSION_STRING << std::endl;
        std::cout << "Build time: " << __TIME__ << " " << __DATE__ << "." << std::endl;
        return false;
    }
    // if no args, print usage.
    std::cerr << parser;
    return false;
}

void MISAMD::onCreate() {
    wall_clock_begin = MPI_Wtime(); // init time recorder.
    ConfigParser *pConfig;
    if (kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) {
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
#ifdef MD_DEV_MODE
// print configure.
    if (kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) {
        std::cout << pConfig->configValues;
    }
#endif

    // prepare logs.
    if (pConfig->configValues.output.logs_mode == LOGS_MODE_CONSOLE && istty()) {
        // set colorful log if we output to console and it is a real tty(no io redirection).
        kiwi::logs::setColorFul(true);
    } else if (pConfig->configValues.output.logs_mode == LOGS_MODE_FILE) {
        kiwi::logs::setLogFile(pConfig->configValues.output.logs_filename);
    }

    // set simulation domain
    MPIDomain::sim_processor = kiwi::mpiUtils::global_process;
    archEnvInit(); // initialize architectures environment.
}

bool MISAMD::prepare() {
    kiwi::logs::d(MASTER_PROCESSOR, "domain", "ranks {}\n", MPIDomain::sim_processor.all_ranks);

    mpi_types::setInterMPIType();
    pSimulation = new MDSimulation(&ConfigParser::getInstance()->configValues);
    const ConfigValues config = ConfigParser::getInstance()->configValues;
    pSimulation->createDomain(config.phaseSpace, config.latticeConst, config.cutoffRadiusFactor); // 区域分解
    // todo alloy ratio seed is not used.
    pSimulation->createAtoms(config.phaseSpace, config.latticeConst, config.timeStepLength,
                             config.createPhaseMode, config.createTSet, config.createSeed, config.alloyRatio);
    return true;
}

void MISAMD::onStart() {
    const ConfigValues config = ConfigParser::getInstance()->configValues;
    pSimulation->prepareForStart(config.potentialFilename);
    kiwi::logs::v(MASTER_PROCESSOR, "simulation", "Start simulation.\n");
    pSimulation->simulate(config.timeSteps); // start simulation.
}

void MISAMD::onFinish() {
    kiwi::logs::s(MASTER_PROCESSOR, "simulation", "finalizing simulation\n");
    //模拟结束
    pSimulation->finalize();
    mpi_types::unsetInterMPIType();

    double wall_clock_end = MPI_Wtime();
    kiwi::logs::i(MASTER_PROCESSOR, "simulation", "total wall clock of program: {}\n",
                  wall_clock_end - wall_clock_begin);
}

void MISAMD::beforeDestroy() {
    kiwi::logs::v(MASTER_PROCESSOR, "app", "app was detached.\n");
    archEnvFinalize(); // clean architectures environment.
}

void MISAMD::onDestroy() {
    delete pSimulation;
}
