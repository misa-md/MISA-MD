#include <iostream>
#include <fstream>
#include "config.h"
#include "utils/mpi_utils.h"

//
// Created by gensh(genshenchu@gmail.com) on 2017/4/16.
//

using namespace std;
config *config::m_pInstance = nullptr;

// a simple single mode.
config *config::newInstance() {
    if (m_pInstance == nullptr) {
        m_pInstance = new config();
    }
    return m_pInstance; // make sure there is a configure instance.
}

config *config::newInstance(const string &configureFilePath) {
    if (m_pInstance == nullptr) {
        m_pInstance = new config(configureFilePath);
    }
    return m_pInstance;
}

config::config() : hasError(false) {
    configValues = new ConfigValues(1024); // todo cap
}

config::config(const string &configurePath) : config() {
    resolveConfig(configurePath);
}

bool config::configureCheck() {
    //todo
    return true;
}

void config::resolveConfig(const string &configurePath) {
// Parse foo.toml. If foo.toml is valid, pr.valid() should be true.
// If not valid, pr.errorReason will contain the parser error reason.
    std::ifstream ifs(configurePath);
    if (!ifs.good()) {
        errorMessage = "can not access the configure file";
        hasError = true;
        return;
    }
    toml::ParseResult pr = toml::parse(ifs);

    if (!pr.valid()) {
        errorMessage = pr.errorReason;
        hasError = true;
        return;
    }
    const toml::Value &v = pr.value;

    // simulation section of config file.
    if (!hasError) {
        resolveConfigSimulation(v);
    }

    // output section of config file.
    if (!hasError) {
        resolveConfigOutput(v);
    }
}

void config::resolveConfigSimulation(const toml::Value &v) {
    // resolve simulation.phasespace
    const toml::Value *tomlPhaseSpace = v.find("simulation.phasespace");
    if (tomlPhaseSpace && tomlPhaseSpace->is<toml::Array>()) {
        const toml::Array &ar = tomlPhaseSpace->as<toml::Array>();
        int index = 0;
        for (const toml::Value &vPS : ar) {
            if (index < DIMENSION && vPS.is<int>()) { //the array index must be less than or equal 3
                configValues->phaseSpace[index] = vPS.as<int>();
            }
            index++;
        }
        if (index != DIMENSION) { //the array length must be 3.
            errorMessage = "array length of value \"simulation.phasespace\" must be 3.";
            hasError = true;
            return;
        }
    }

    //resolve simulation.cutoff_radius
    const toml::Value *tomlCutoffRadius = v.find("simulation.cutoff_radius");
    if (tomlCutoffRadius && tomlCutoffRadius->is<double>()) {
        configValues->cutoffRadius = tomlCutoffRadius->as<double>();
    }

    //resolve simulation.latticeconst
    const toml::Value *tomlLatticeConst = v.find("simulation.lattice_const");
    if (tomlLatticeConst && tomlLatticeConst->is<double>()) {
        configValues->latticeConst = tomlLatticeConst->as<double>();
    }

    //resolve simulation.timesteps
    const toml::Value *tomlTimeSteps = v.find("simulation.timesteps");
    if (tomlTimeSteps && tomlTimeSteps->is<long>()) {
        configValues->timeSteps = tomlTimeSteps->as<long>();
    }

    //resolve simulation.createphase
    const toml::Value *tomlCreatePhase = v.find("simulation.createphase.create_phase");
    if (tomlCreatePhase && tomlCreatePhase->is<bool>()) {
        configValues->createPhaseMode = tomlCreatePhase->as<bool>();
        if (configValues->createPhaseMode) { //create mode
            const toml::Value *tomlTSet = v.find("simulation.createphase.create_t_set");
            if (tomlTSet && tomlTSet->is<double>()) {
                configValues->createTSet = tomlTSet->as<double>();
            }
            const toml::Value *tomlSeed = v.find("simulation.createphase.create_seed");
            if (tomlSeed && tomlSeed->is<int>()) {
                configValues->createSeed = tomlSeed->as<int>();
            }
        } else {  //read mode.
            const toml::Value *tomlTSet = v.find("simulation.createphase.read_phase_filename");
            if (tomlTSet && tomlTSet->is<string>()) {
                configValues->readPhaseFilename = tomlTSet->as<string>();
            } else {
                errorMessage = "read phase file must be specified.";
                hasError = true;
                return;
            }
        }
    } else {
        errorMessage = "create phase mode(read/create) is required.";
        hasError = true;
        return;
    }

    //resolve simulation.collision
    const toml::Value *tomlCollisionSteps = v.find("simulation.collision.collision_steps");
    if (tomlCollisionSteps && tomlCollisionSteps->is<long>()) {
        configValues->collisionSteps = tomlCollisionSteps->as<long>();
    }
    const toml::Value *tomlCollisionLat = v.find("simulation.collision.lat");
    if (tomlCollisionLat && tomlCollisionLat->is<toml::Array>()) {
        const toml::Array &ar = tomlCollisionLat->as<toml::Array>();
        int index = 0;
        for (const toml::Value &value : ar) {
            if (index < 4 && value.is<int>()) { //the array index must be less than or equal 4
                configValues->collisionLat[index] = value.as<int>();
            }
            index++;
        }
        if (index != 4) { //the array length must be 3.
            errorMessage = "array length of value \"simulation..collision.lat\" must be 4.";
            hasError = true;
            return;
        }
    }
    const toml::Value *tomlCollisionV = v.find("simulation.collision.collision_v");
    if (tomlCollisionV && tomlCollisionV->is<toml::Array>()) {
        const toml::Array &ar = tomlCollisionV->as<toml::Array>();
        int index = 0;
        for (const toml::Value &value : ar) {
            if (index < DIMENSION && value.is<double>()) { //the array index must be less than or equal 3
                configValues->collisionV[index] = value.as<double>();
            }
            index++;
        }
        if (index != DIMENSION) { //the array length must be 3.
            errorMessage = "array length of value \"simulation..collision.collision_v\" must be 3.";
            hasError = true;
            return;
        }
    }

    //potential_file
    const toml::Value *tomlPotentialFileType = v.find("simulation.potential_file.type");
    if (tomlPotentialFileType && tomlPotentialFileType->is<string>()) {
        configValues->potentialFileType = tomlPotentialFileType->as<string>();
    } else {
        errorMessage = "potential file type must be specified.";
        hasError = true;
        return;
    }
    const toml::Value *tomlPotentialFilename = v.find("simulation.potential_file.filename");
    if (tomlPotentialFilename && tomlPotentialFilename->is<string>()) {
        configValues->potentialFilename = tomlPotentialFilename->as<string>();
    } else {
        errorMessage = "potential filename must be specified.";
        hasError = true;
        return;
    }
}

void config::resolveConfigOutput(const toml::Value &v) {
    const toml::Value *tomlOutputMode = v.find("output.mode");
    if (tomlOutputMode && tomlOutputMode->is<string>()) {
        if (tomlOutputMode->as<string>() == "copy") {
            configValues->outputMode = OUTPUT_COPY_MODE;
        } else {
            configValues->outputMode = OUTPUT_DIRECT_MODE;
        }
    }

    // todo check if it is a path.
    const toml::Value *tomlOutputDumpFilename = v.find("output.dump_filename");
    if (tomlOutputDumpFilename && tomlOutputDumpFilename->is<string>()) {
        configValues->outputDumpFilename = tomlOutputDumpFilename->as<string>();
    } else {
        configValues->outputDumpFilename = DEFAULT_OUTPUT_DUMP_FILENAME;
    }
}

void config::sync() {
    configValues->newPackBuffer();
    if (mpiUtils::ownRank == MASTER_PROCESSOR) { // pack data.
        configValues->packdata();
    }

    MPI_Bcast(configValues->getPackedData(), configValues->getPackedDataCap(),
              MPI_BYTE, MASTER_PROCESSOR, MPI_COMM_WORLD); // synchronize config information

    if (mpiUtils::ownRank != MASTER_PROCESSOR) { // unpack data.
        configValues->unpackdata();
    }
    configValues->releasePackBuffer(); // release memory after usage.
}
