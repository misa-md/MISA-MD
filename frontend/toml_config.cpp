//
// Created by gensh(genshenchu@gmail.com) on 2017/4/16.
//

#include <iostream>
#include <fstream>
#include <utils/mpi_utils.h>
#include <utils/bundle.h>
#include <climits>

#include "toml_config.h"
#include "utils/rpcc.hpp"
#include "def_config_values.h"


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
    configValues.timeStepLength = table->get_as<double>("def_timesteps_length").value_or(default_time_length);
    resolveVariableStepLen(table);

    //resolve simulation.createphase
    auto tableCreatephase = table->get_table("createphase");

    configValues.createPhaseMode = tableCreatephase->get_as<bool>("create_phase")
            .value_or(default_create_phase); // default value is true
    if (configValues.createPhaseMode) { //create mode
        auto tomlTSet = tableCreatephase->get_as<double>("create_t_set");
        if (tomlTSet) {
            configValues.createTSet = *tomlTSet;
        } else {
            setError("creation t_set must be specified.");
            return;
        }
        configValues.createSeed = tableCreatephase->get_as<int>("create_seed").value_or(default_random_seek);
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

void ConfigParser::resolveVariableStepLen(std::shared_ptr<cpptoml::table> table) {
    auto nested = table->get_array_of<cpptoml::array>("variable_step_length");
    if (!nested) { // if it is not set.
        return;
    }
    cpptoml::option<std::vector<int64_t>> break_steps_i64 = (*nested)[0]->get_array_of<int64_t>();
    cpptoml::option<std::vector<double>> step_lens = (*nested)[1]->get_array_of<double>();
    if (break_steps_i64->size() != step_lens->size()) {
        setError("size of break points and step length is not equal in variable timestep length configuration.");
        return;
    }
    const unsigned long _size = break_steps_i64->size();
    std::vector<unsigned long> break_steps;
    configValues.vsl_size = _size; // set size
    for (const int64_t &val : *break_steps_i64) { // type conversion from int64 to unsigned long.
        break_steps.push_back(val);
    }
    configValues.setVarStepLengths(break_steps, *step_lens, _size); // set variable length array.
}

void ConfigParser::resolveConfigOutput(std::shared_ptr<cpptoml::table> table) {
    auto tomlAtomsDumpMode = table->get_as<std::string>("atoms_dump_mode");
    if (tomlAtomsDumpMode) {
        if ("copy" == *tomlAtomsDumpMode) { // todo equal?
            configValues.output.atomsDumpMode = OutputMode::COPY;
        } else {
            configValues.output.atomsDumpMode = OutputMode::DEBUG;
        }
    }
    configValues.output.atomsDumpInterval = table->get_as<uint64_t>("atoms_dump_interval").value_or(ULONG_MAX);
    configValues.output.outByFrame = table->get_as<bool>("by_frame").value_or(false);

    configValues.output.atomsDumpFilePath = table->get_as<std::string>("atoms_dump_file_path")
            .value_or(DEFAULT_OUTPUT_DUMP_FILE_PATH);
    configValues.output.originDumpPath = table->get_as<std::string>("origin_dump_path").value_or("");
    // todo check if it is a real path.
    if (configValues.output.outByFrame && configValues.output.atomsDumpFilePath.find("{}") == std::string::npos) {
        setError("error format of dump file path");
    }

    // resolve logs.
    std::shared_ptr<cpptoml::table> logs_table = table->get_table("logs");
    auto logsModeString = logs_table->get_as<std::string>("logs_mode").value_or(DEFAULT_LOGS_MODE_CONSOLE_STRING);
    if (LOGS_MODE_CONSOLE_STRING == logsModeString) {
        configValues.output.logs_mode = LOGS_MODE_CONSOLE;
    } else {
        configValues.output.logs_mode = LOGS_MODE_FILE;
        configValues.output.logs_filename = logs_table->get_as<std::string>("logs_filename").value_or("");
        // if filename is empty, we will generate a filename.
        if (configValues.output.logs_filename.empty()) {
            std::ostringstream str_stream;
            str_stream << "md.";
            str_stream << rpcc();
            str_stream << ".log";
            configValues.output.logs_filename = str_stream.str();
        }
    }
}

void ConfigParser::resolveConfigAlloy(std::shared_ptr<cpptoml::table> table) {
    configValues.alloyCreateSeed = table->get_as<int>("create_seed").value_or(default_random_seek);
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
