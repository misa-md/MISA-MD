//
// Created by genshen(genshenchu@gmail.com) on 2018-3-9.
//

#include <iostream>
#include <cassert>
#include <mpi.h>
#include <utils/bundle.h>
#include "config_values.h"
#include "def_config_values.h"

Stage::Stage() : steps(0), step_length(default_time_length),
                 collisionStep(0), collisionLat{0, 0, 0, 0},
                 pkaEnergy(0), direction{1.0, 1.0, 1.0} {}

void Stage::packdata(kiwi::Bundle &bundle) {
    bundle.put(steps);
    bundle.put(step_length);
    // collision
    bundle.put(collision_set);
    bundle.put(collisionStep);
    bundle.put(4, collisionLat);
    bundle.put(pkaEnergy);
    bundle.put(DIMENSION, direction);
}

void Stage::unnpackdata(int &cursor, kiwi::Bundle &bundle) {
    bundle.get(cursor, steps); // get size,
    bundle.get(cursor, step_length);
    // collision
    bundle.get(cursor, collision_set);
    bundle.get(cursor, collisionStep);
    bundle.get(cursor, 4, collisionLat);
    bundle.get(cursor, pkaEnergy);
    bundle.get(cursor, DIMENSION, direction);
}

ConfigValues::ConfigValues() :
        phaseSpace{0, 0, 0}, cutoffRadiusFactor(0.0), latticeConst(0.0),
        timeSteps(0),
        createPhaseMode(true), createTSet(0.0), createSeed(1024), readPhaseFilename(""),
        alloyCreateSeed(1024), alloyRatio{1, 0, 0},
        output(), stages() {}
// todo potential type and filename initialize.

void ConfigValues::packdata(kiwi::Bundle &bundle) {
    // append data into buffer.
    bundle.put(DIMENSION, phaseSpace); // todo remove MPI_COMM_WORLD to initial method.
    bundle.put(cutoffRadiusFactor);
    bundle.put(latticeConst);

    // step and step length
    bundle.put(timeSteps);
    bundle.put(timeStepLength);

    bundle.put(createPhaseMode);
    bundle.put(createTSet);
    bundle.put(createSeed);
    bundle.put(readPhaseFilename);

    // alloy
    bundle.put(alloyCreateSeed);
    bundle.put(atom_type::num_atom_types, alloyRatio);

    bundle.put(potentialFileType);
    bundle.put(potentialFilename);

    // output section
    bundle.put(output.atomsDumpMode);
    bundle.put(output.atomsDumpFilePath);
    bundle.put(output.originDumpPath);
    bundle.put(output.atomsDumpInterval);
    bundle.put(output.outByFrame);

    bundle.put(output.thermo_interval);
    // logs subsection in output section.
    bundle.put(output.logs_mode);
    bundle.put(output.logs_filename);

    // stages
    bundle.put(stages.size());
    for (Stage stage:stages) {
        stage.packdata(bundle);
    }
}

void ConfigValues::unpackdata(kiwi::Bundle &bundle) {
    // fetch data from buffer.
//    if (getPackedData() != nullptr) { // buffer != null
    int cursor = 0;

    bundle.get(cursor, DIMENSION, phaseSpace);
    bundle.get(cursor, cutoffRadiusFactor);
    bundle.get(cursor, latticeConst);

    // step and step length
    bundle.get(cursor, timeSteps);
    bundle.get(cursor, timeStepLength);

    bundle.get(cursor, createPhaseMode);
    bundle.get(cursor, createTSet);
    bundle.get(cursor, createSeed);
    bundle.get(cursor, readPhaseFilename);

    // alloy
    bundle.get(cursor, alloyCreateSeed);
    bundle.get(cursor, atom_type::num_atom_types, alloyRatio);

    bundle.get(cursor, potentialFileType);
    bundle.get(cursor, potentialFilename);

    // output section.
    bundle.get(cursor, output.atomsDumpMode);
    bundle.get(cursor, output.atomsDumpFilePath);
    bundle.get(cursor, output.originDumpPath);
    bundle.get(cursor, output.atomsDumpInterval);
    bundle.get(cursor, output.outByFrame);

    bundle.get(cursor, output.thermo_interval);

    // logs subsection in output section.
    bundle.get(cursor, output.logs_mode);
    bundle.get(cursor, output.logs_filename);

    // stages
    std::size_t stages_size = 0;
    bundle.get(cursor, stages_size);
    stages.resize(stages_size);
    for (std::size_t i = 0; i < stages_size; i++) {
        stages[i].unnpackdata(cursor, bundle);
    }
}

std::ostream &operator<<(std::ostream &os, const ConfigValues &cv) {
    // simulation section
    os << "===========config of simulation=============" << std::endl;
    os << "simulation.phase_space:" << cv.phaseSpace[0] << "," << cv.phaseSpace[1]
       << "," << cv.phaseSpace[2] << "," << std::endl;
    os << "simulation.cutoff_radius:" << cv.cutoffRadiusFactor << std::endl;
    os << "simulation.lattice_const:" << cv.latticeConst << std::endl;
    os << "simulation.timesteps:" << cv.timeSteps << std::endl;

    os << "simulation.createphase.createPhaseMode:" << (cv.createPhaseMode ? "true" : "false") << std::endl;
    os << "simulation.createphase.tSet:" << cv.createTSet << std::endl;
    os << "simulation.createphase.seed:" << cv.createSeed << std::endl;
    os << "simulation.createphase.readPhaseFilename:" << cv.readPhaseFilename << std::endl;

    // simulation.alloy
    os << "simulation.alloy.alloyCreateSeed:" << cv.alloyCreateSeed << std::endl;
    os << "simulation.alloy.ratio Fe:Cu:Ni\t" << cv.alloyRatio[atom_type::Fe] <<
       ":" << cv.alloyRatio[atom_type::Cu] << ":" << cv.alloyRatio[atom_type::Ni] << ":" << std::endl;

    os << "simulation.potential_file.type:" << cv.potentialFileType << std::endl;
    os << "simulation.potential_file.filename:" << cv.potentialFilename << std::endl;

    // output section
    os << "output.mode(debug:0, copy:1):" << cv.output.atomsDumpMode << std::endl;
    os << "output.dump_interval:" << cv.output.atomsDumpInterval << std::endl;
    os << "output.dump_file_path:" << cv.output.atomsDumpFilePath << std::endl;
    os << "output.origin_dump_path:" << cv.output.originDumpPath << std::endl;
    os << "output.thermo.interval:" << cv.output.thermo_interval << std::endl;
    os << "output.logs.mode: "
       << (cv.output.logs_mode == LOGS_MODE_CONSOLE ? LOGS_MODE_CONSOLE_STRING : LOGS_MODE_FILE_STRING)
       << ", output.logs.by-frame:" << cv.output.outByFrame << std::endl;
    os << "output.logs_filename:" << cv.output.logs_filename << std::endl;

    // stages
    os << "stages:" << std::endl;
    for (Stage stage : cv.stages) {
        os << "stage.steps:" << stage.steps << std::endl;
        os << "stage.steps_length:" << stage.step_length << std::endl;
        if (stage.collision_set) {
            os << "stage.collision.collision_step:" << stage.collisionStep << std::endl;
            os << "stage.collision.lat:" << stage.collisionLat[0] << "," << stage.collisionLat[1] << ","
               << stage.collisionLat[2] << "," << stage.collisionLat[3] << std::endl;
            os << "pka:" << stage.pkaEnergy << ",";
            os << "direction:" << stage.direction[0] << "," << stage.direction[1] << ","
               << stage.direction[2] << std::endl;
        }
    }
    os << "============================================" << std::endl << std::endl;
    return os;
}
