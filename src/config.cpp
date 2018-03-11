#include <iostream>
#include <fstream>
#include "config.h"
#include "toml.hpp"
#include "pre_config.h"

//
// Created by gensh(genshenchu@gmail.com) on 2017/4/16.
//

using namespace std;
config *config::m_pInstance = nullptr;

//single mode
config *config::newInstance(string configureFilePath) {
    if (m_pInstance == nullptr) {
        m_pInstance = new config(configureFilePath);
    }
    return m_pInstance;
}

config *config::newInstance() {
    if (m_pInstance == nullptr) {
        return new config();
    }
    return m_pInstance; // make sure there is a configure instance.
}

config::config() : phaseSpace{0, 0, 0}, cutoffRadius(0.0), latticeConst(0.0), timeSteps(10),
                   createPhaseMode(true), createTSet(0.0), createSeed(1), readPhaseFilename(""),
                   collisionSteps(0), collisionLat{0, 0, 0, 0}, collisionV{0.0, 0.0, 0.0},
                   hasError(false) {
}

config::config(string configurePath) : config() {
    resolveConfig(configurePath);
}

void config::onPostMPICopy(config *ptrConfig) {
    m_pInstance = ptrConfig;
}

bool config::configureCheck() {
    //todo
    return true;
}

void config::resolveConfig(string configurePath) {
// Parse foo.toml. If foo.toml is valid, pr.valid() should be true.
// If not valid, pr.errorReason will contain the parser error reason.
    std::ifstream ifs(configurePath);
    toml::ParseResult pr = toml::parse(ifs);

    if (!pr.valid()) {
        errorMessage = pr.errorReason;
        hasError = true;
        return;
    }
    const toml::Value &v = pr.value;

    //resolve simulation.phasespace
    const toml::Value *tomlPhaseSpace = v.find("simulation.phasespace");
    if (tomlPhaseSpace && tomlPhaseSpace->is<toml::Array>()) {
        const toml::Array &ar = tomlPhaseSpace->as<toml::Array>();
        int index = 0;
        for (const toml::Value &vPS : ar) {
            if (index < DIMENSION && vPS.is<int>()) { //the array index must be less than or equal 3
                phaseSpace[index] = vPS.as<int>();
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
        cutoffRadius = tomlCutoffRadius->as<double>();
    }

    //resolve simulation.latticeconst
    const toml::Value *tomlLatticeConst = v.find("simulation.lattice_const");
    if (tomlLatticeConst && tomlLatticeConst->is<double>()) {
        latticeConst = tomlLatticeConst->as<double>();
    }

    //resolve simulation.timesteps
    const toml::Value *tomlTimeSteps = v.find("simulation.timesteps");
    if (tomlTimeSteps && tomlTimeSteps->is<long>()) {
        timeSteps = tomlTimeSteps->as<long>();
    }

    //resolve simulation.createphase
    const toml::Value *tomlCreatePhase = v.find("simulation.createphase.create_phase");
    if (tomlCreatePhase && tomlCreatePhase->is<bool>()) {
        createPhaseMode = tomlCreatePhase->as<bool>();
        if (createPhaseMode) { //create mode
            const toml::Value *tomlTSet = v.find("simulation.createphase.create_t_set");
            if (tomlTSet && tomlTSet->is<double>()) {
                createTSet = tomlTSet->as<double>();
            }
            const toml::Value *tomlSeed = v.find("simulation.createphase.create_seed");
            if (tomlSeed && tomlSeed->is<int>()) {
                createSeed = tomlSeed->as<int>();
            }
        } else {  //read mode.
            const toml::Value *tomlTSet = v.find("simulation.createphase.read_phase_filename");
            if (tomlTSet && tomlTSet->is<string>()) {
                readPhaseFilename = tomlTSet->as<string>();
            } else {
                errorMessage = "read phase file must be specified.";
                hasError = true;
                return;
            }
        }
    } else {
        errorMessage = "create phase mode(read/create) is required..";
        hasError = true;
        return;
    }

    //resolve simulation.collision
    const toml::Value *tomlCollisionSteps = v.find("simulation.collision.collision_steps");
    if (tomlCollisionSteps && tomlCollisionSteps->is<long>()) {
        collisionSteps = tomlCollisionSteps->as<long>();
    }
    const toml::Value *tomlCollisionLat = v.find("simulation.collision.lat");
    if (tomlCollisionLat && tomlCollisionLat->is<toml::Array>()) {
        const toml::Array &ar = tomlCollisionLat->as<toml::Array>();
        int index = 0;
        for (const toml::Value &value : ar) {
            if (index < 4 && value.is<int>()) { //the array index must be less than or equal 4
                collisionLat[index] = value.as<int>();
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
                collisionV[index] = value.as<double>();
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
        potentialFileType = tomlPotentialFileType->as<string>();
    } else {
        errorMessage = "potential file type must be specified.";
        hasError = true;
        return;
    }
    const toml::Value *tomlPotentialFilename = v.find("simulation.potential_file.filename");
    if (tomlPotentialFilename && tomlPotentialFilename->is<string>()) {
        potentialFilename = tomlPotentialFilename->as<string>();
    } else {
        errorMessage = "potential filename must be specified.";
        hasError = true;
        return;
    }

#ifdef DEV_MODE
    cout << "config of simulation:" << endl;
    cout << "simulation.phase_space:" << phaseSpace[0] << phaseSpace[1] << phaseSpace[2] << endl;
    cout << "simulation.cutoff_radius:" << cutoffRadius << endl;
    cout << "simulation.lattice_const:" << latticeConst << endl;
    cout << "simulation.timesteps:" << timeSteps << endl;

    cout << "simulation.createphase.tSet:" << createTSet << endl;
    cout << "simulation.createphase.seed:" << createSeed << endl;

    cout << "simulation.collision.collision_steps:" << collisionSteps << endl;
    cout << "simulation.collision.lat:" << collisionLat[0] << collisionLat[1] << collisionLat[2] << collisionLat[3]
         << endl;
    cout << "simulation.collision.collision_v:" << collisionV[0] << collisionV[1] << collisionV[2] << endl;

    cout << "simulation.potential_file.type:" << potentialFileType << endl;
    cout << "simulation.potential_file.filename:" << potentialFilename << endl;

#endif
}
