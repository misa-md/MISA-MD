#include <iostream>
#include <fstream>
#include <utils/mpi_utils.h>
#include <utils/bundle.h>

#include "toml_config.h"

//
// Created by gensh(genshenchu@gmail.com) on 2017/4/16.
//

using namespace std;
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

ConfigParser *ConfigParser::newInstance(const string &configureFilePath) {
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

void ConfigParser::resolveConfigSimulation(std::shared_ptr<cpptoml::table> v) {
    // resolve simulation.phasespace
    auto tomlPhaseSpace = v->get_array_of<int64_t>("phasespace");
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
    auto tomlCutoffRadiusFactor = v->get_as<double>("cutoff_radius_factor");
    if (tomlCutoffRadiusFactor) {
        configValues.cutoffRadiusFactor = *tomlCutoffRadiusFactor;
    }

    //resolve simulation.latticeconst
    auto tomlLatticeConst = v->get_as<double>("lattice_const");
    if (tomlLatticeConst) {
        configValues.latticeConst = *tomlLatticeConst;
    }

    //resolve simulation.timesteps
    auto tomlTimeSteps = v->get_as<unsigned long>("timesteps");
    if (tomlTimeSteps) {
        configValues.timeSteps = *tomlTimeSteps;
    }

    //resolve simulation.createphase
    auto tableCreatephase = v->get_table("createphase");

    auto tomlCreatePhase = tableCreatephase->get_as<bool>("create_phase");
    if (tomlCreatePhase) { // todo bool.
        configValues.createPhaseMode = *tomlCreatePhase;

        if (configValues.createPhaseMode) { //create mode
            auto tomlTSet = tableCreatephase->get_as<double>("create_t_set");
            if (tomlTSet) {
                configValues.createTSet = *tomlTSet;
            }
            auto tomlSeed = tableCreatephase->get_as<int>("create_seed");
            if (tomlSeed) {
                configValues.createSeed = *tomlSeed;
            }
        } else {  //read mode.
            auto tomlTSet = tableCreatephase->get_as<std::string>("read_phase_filename");
            if (tomlTSet) {
                configValues.readPhaseFilename = *tomlTSet;
            } else {
                setError("read phase file must be specified.");
                return;
            }
        }
    } else {
        setError("create phase mode(read/create) is required.");
        return;
    }

    //resolve simulation.collision
    auto tomlCollision = v->get_table("collision");

    auto tomlCollisionSteps = tomlCollision->get_as<unsigned long>("collision_steps");
    if (tomlCollisionSteps) {
        configValues.collisionSteps = *tomlCollisionSteps;
    }
    auto tomlCollisionLat = tomlCollision->get_array_of<int64_t>("lat");
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
            setError("array length of value \"simulation..collision.lat\" must be 4.");
            return;
        }
    }
    auto tomlCollisionV = tomlCollision->get_array_of<double>("collision_v");
    if (tomlCollisionV) {
        int index = 0;
        for (auto &value : *tomlCollisionV) {
            if (index < DIMENSION) { //the array index must be less than or equal 3
                configValues.collisionV[index] = value;
            }
            index++;
        }
        if (index != DIMENSION) { //the array length must be 3.
            setError("array length of value \"simulation..collision.collision_v\" must be 3.");
            return;
        }
    }

    //potential_file
    auto tomlPotentialFile = v->get_table("potential_file");
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

void ConfigParser::resolveConfigOutput(shared_ptr<cpptoml::table> v) {
    auto tomlOutputMode = v->get_as<std::string>("mode");
    if (tomlOutputMode) {
        if ("copy" == *tomlOutputMode) { // todo equal?
            configValues.outputMode = OUTPUT_COPY_MODE;
        } else {
            configValues.outputMode = OUTPUT_DIRECT_MODE;
        }
    }

    // todo check if it is a path.
    auto tomlOutputDumpFilename = v->get_as<std::string>("dump_filename");
    if (tomlOutputDumpFilename) {
        configValues.outputDumpFilename = *tomlOutputDumpFilename;
    } else {
        configValues.outputDumpFilename = DEFAULT_OUTPUT_DUMP_FILENAME;
    }
}
