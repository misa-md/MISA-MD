#include <iostream>
#include <fstream>
#include <utils/mpi_utils.h>
#include <utils/bundle.h>
#include <climits>

#include "toml_config.h"
#include "utils/rpcc.hpp"

//
// Created by gensh(genshenchu@gmail.com) on 2017/4/16.
//

ConfigParser *ConfigParser::m_pInstance = nullptr;

ConfigParser::ConfigParser() : kiwi::config::config() {
    configValues = ConfigValues(); // todo ConfigValues cap
}

// a simple single mode.
ConfigParser *ConfigParser::getInstance() {
    if (m_pInstance == nullptr) {
        m_pInstance = new ConfigParser();
    }
    return m_pInstance; // make sure there is a configure instance.
}

ConfigParser *ConfigParser::newInstance(const std::string &configureFilePath) {
    if (m_pInstance == nullptr) {
        m_pInstance = new ConfigParser();  // todo delete
        m_pInstance->resolve(configureFilePath);
    }
    return m_pInstance;
}

bool ConfigParser::configureCheck() {
    //todo
    return true;
}

// @override only for master processor.
// @see https://github.com/skystrife/cpptoml#example-usage for more details.
void ConfigParser::resolveConfig(std::shared_ptr<cpptoml::table> table) {
    resolveConfigSimulation(table->get_table("simulation"));
    resolveConfigOutput(table->get_table("output"));
}

// @override
void ConfigParser::putConfigData(kiwi::Bundle &bundle) {
    configValues.packdata(bundle);
}

// @override
void ConfigParser::getConfigData(kiwi::Bundle &bundle) {
    configValues.unpackdata(bundle);
}

void ConfigParser::resolveConfigSimulation(std::shared_ptr<cpptoml::table> table) {
    // resolve simulation.phasespace
    auto tomlPhaseSpace = table->get_array_of<int64_t>("phasespace");
    if (tomlPhaseSpace) {
        int index = 0;
        for (const auto &val : *tomlPhaseSpace) {
            if (index < DIMENSION) { //the array index must be less than or equal 3
                configValues.phaseSpace[index] = val;  // todo int64_t
            }
            index++;
        }
        if (index != DIMENSION) { //the array length must be 3.
            setError("array length of value \"simulation.phasespace\" must be 3.");
            return;
        }
    }

    //resolve simulation.cutoff_radius
    auto tomlCutoffRadiusFactor = table->get_as<double>("cutoff_radius_factor");
    if (tomlCutoffRadiusFactor) {
        configValues.cutoffRadiusFactor = *tomlCutoffRadiusFactor;
    }

    //resolve simulation.latticeconst
    auto tomlLatticeConst = table->get_as<double>("lattice_const");
    if (tomlLatticeConst) {
        configValues.latticeConst = *tomlLatticeConst;
    }

    //resolve simulation.timesteps
    auto tomlTimeSteps = table->get_as<unsigned long>("timesteps");
    if (tomlTimeSteps) {
        configValues.timeSteps = *tomlTimeSteps;
    }
    configValues.timeStepLength = table->get_as<double>("timesteps_length").value_or(const_default_time_length);

    //resolve simulation.createphase
    auto tableCreatephase = table->get_table("createphase");

    configValues.createPhaseMode = tableCreatephase->get_as<bool>("create_phase")
            .value_or(const_default_create_phase); // default value is true
    if (configValues.createPhaseMode) { //create mode
        auto tomlTSet = tableCreatephase->get_as<double>("create_t_set");
        if (tomlTSet) {
            configValues.createTSet = *tomlTSet;
        } else {
            setError("creation t_set must be specified.");
            return;
        }
        configValues.createSeed = tableCreatephase->get_as<int>("create_seed").value_or(const_default_random_seek);
    } else {  // read mode.
        auto tomlTSet = tableCreatephase->get_as<std::string>("read_phase_filename");
        if (tomlTSet) {
            configValues.readPhaseFilename = *tomlTSet;
        } else {
            setError("read phase file must be specified.");
            return;
        }
    }

    // resolve simulation.alloy
    resolveConfigAlloy(table->get_table("alloy"));
    //resolve simulation.collision
    resolveConfigCollision(table->get_table("collision"));

    //potential_file
    auto tomlPotentialFile = table->get_table("potential_file");
    auto tomlPotentialFileType = tomlPotentialFile->get_as<std::string>("type");
    if (tomlPotentialFileType) {
        configValues.potentialFileType = *tomlPotentialFileType;
    } else {
        setError("potential file type must be specified.");
        return;
    }
    auto tomlPotentialFilename = tomlPotentialFile->get_as<std::string>("filename");
    if (tomlPotentialFilename) {
        configValues.potentialFilename = *tomlPotentialFilename;
    } else {
        setError("potential filename must be specified.");
        return;
    }
}

void ConfigParser::resolveConfigOutput(std::shared_ptr<cpptoml::table> table) {
    auto tomlAtomsDumpMode = table->get_as<std::string>("atoms_dump_mode");
    if (tomlAtomsDumpMode) {
        if ("copy" == *tomlAtomsDumpMode) { // todo equal?
            configValues.atomsDumpMode = OUTPUT_COPY_MODE;
        } else {
            configValues.atomsDumpMode = OUTPUT_DIRECT_MODE;
        }
    }
    configValues.atomsDumpInterval = table->get_as<uint64_t>("atoms_dump_interval").value_or(ULONG_MAX);
    configValues.outByFrame = table->get_as<bool>("by_frame").value_or(false);

    configValues.atomsDumpFilePath = table->get_as<std::string>("atoms_dump_file_path")
            .value_or(DEFAULT_OUTPUT_DUMP_FILE_PATH);
    configValues.originDumpPath = table->get_as<std::string>("origin_dump_path").value_or("");
    // todo check if it is a real path.
    if (configValues.outByFrame && configValues.atomsDumpFilePath.find("{}") == std::string::npos) {
        setError("error format of dump file path");
    }

    // resolve logs.
    std::shared_ptr<cpptoml::table> logs_table = table->get_table("logs");
    auto logsModeString = logs_table->get_as<std::string>("logs_mode").value_or(DEFAULT_LOGS_MODE_CONSOLE_STRING);
    if (LOGS_MODE_CONSOLE_STRING == logsModeString) {
        configValues.logs_mode = LOGS_MODE_CONSOLE;
    } else {
        configValues.logs_mode = LOGS_MODE_FILE;
        configValues.logs_filename = logs_table->get_as<std::string>("logs_filename").value_or("");
        // if filename is empty, we will generate a filename.
        if (configValues.logs_filename.empty()) {
            std::ostringstream str_stream;
            str_stream << "md.";
            str_stream << rpcc();
            str_stream << ".log";
            configValues.logs_filename = str_stream.str();
        }
    }
}

void ConfigParser::resolveConfigAlloy(std::shared_ptr<cpptoml::table> table) {
    configValues.alloyCreateSeed = table->get_as<int>("create_seed").value_or(const_default_random_seek);
    auto tomlAlloyRatioFe = table->get_qualified_as<int>("ratio.Fe");
    if (tomlAlloyRatioFe) {
        configValues.alloyRatio[atom_type::Fe] = *tomlAlloyRatioFe;
    }
    auto tomlAlloyRatioCu = table->get_qualified_as<int>("ratio.Cu");
    if (tomlAlloyRatioCu) {
        configValues.alloyRatio[atom_type::Cu] = *tomlAlloyRatioCu;
    }
    auto tomlAlloyRatioNi = table->get_qualified_as<int>("ratio.Ni");
    if (tomlAlloyRatioNi) {
        configValues.alloyRatio[atom_type::Ni] = *tomlAlloyRatioNi;
    }
}

void ConfigParser::resolveConfigCollision(std::shared_ptr<cpptoml::table> table) {
    auto tomlCollisionStep = table->get_as<unsigned long>("collision_step");
    if (tomlCollisionStep) {
        configValues.collisionStep = *tomlCollisionStep;
    }
    auto tomlCollisionLat = table->get_array_of<int64_t>("lat");
    if (tomlCollisionLat) {
//        const toml::Array &ar = tomlCollisionLat->as<toml::Array>();
        int index = 0;
        for (auto &value: *tomlCollisionLat) {
            if (index < 4) { //the array index must be less than or equal 4
                configValues.collisionLat[index] = value;  // todo conversion.
            }
            index++;
        }
        if (index != 4) { //the array length must be 3.
            setError("array length of value \"simulation.collision.lat\" must be 4.");
            return;
        }
    }
    configValues.pkaEnergy = table->get_as<double>("pka").value_or(0.0);

    auto tomlCollisionV = table->get_array_of<double>("direction");
    if (tomlCollisionV) {
        int index = 0;
        for (auto &value : *tomlCollisionV) {
            if (index < DIMENSION) { //the array index must be less than or equal 3
                configValues.direction[index] = value;
            }
            index++;
        }
        if (index != DIMENSION) { //the array length must be 3.
            setError("array length of value \"simulation..collision.collision_v\" must be 3.");
            return;
        }
    }
}
