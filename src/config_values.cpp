//
// Created by genshen(genshenchu@gmail.com) on 2018-3-9.
//

#include <mpi.h>
#include <iostream>
#include <utils/bundle.h>
#include "config_values.h"

ConfigValues::ConfigValues() :
        phaseSpace{0, 0, 0}, cutoffRadiusFactor(0.0), latticeConst(0.0), timeSteps(10),
        createPhaseMode(true), createTSet(0.0), createSeed(1024), readPhaseFilename(""),
        alloyCreateSeed(1024), alloyRatio{1, 0, 0},
        collisionStep(0), collisionLat{0, 0, 0, 0}, pkaEnergy(0), direction{1.0, 1.0, 1.0},
        // todo potential type and filename initialize.
        atomsDumpMode(OUTPUT_COPY_MODE), atomsDumpFilePath(DEFAULT_OUTPUT_DUMP_FILE_PATH),
        logs_mode(LOGS_MODE_CONSOLE), logs_filename("") {}

void ConfigValues::packdata(kiwi::Bundle &bundle) {
    // append data into buffer.
    bundle.put(DIMENSION, phaseSpace); // todo remove MPI_COMM_WORLD to initial method.
    bundle.put(cutoffRadiusFactor);
    bundle.put(latticeConst);
    bundle.put(timeSteps);
    bundle.put(timeStepLength);

    bundle.put(createPhaseMode);
    bundle.put(createTSet);
    bundle.put(createSeed);
    bundle.put(readPhaseFilename);

    // alloy
    bundle.put(alloyCreateSeed);
    bundle.put(atom_type::num_atom_types, alloyRatio);

    bundle.put(collisionStep);
    bundle.put(4, collisionLat);
    bundle.put(pkaEnergy);
    bundle.put(DIMENSION, direction);

    bundle.put(potentialFileType);
    bundle.put(potentialFilename);

    // out section
    bundle.put(atomsDumpMode);
    bundle.put(atomsDumpFilePath);
    bundle.put(atomsDumpInterval);
    bundle.put(outByFrame);
    // logs subsection in output section.
    bundle.put(logs_mode);
    bundle.put(logs_filename);
}

void ConfigValues::unpackdata(kiwi::Bundle &bundle) {
    // fetch data from buffer.
//    if (getPackedData() != nullptr) { // buffer != null
    int cursor = 0;
    bundle.get(cursor, DIMENSION, phaseSpace);
    bundle.get(cursor, cutoffRadiusFactor);
    bundle.get(cursor, latticeConst);
    bundle.get(cursor, timeSteps);
    bundle.get(cursor, timeStepLength);

    bundle.get(cursor, createPhaseMode);
    bundle.get(cursor, createTSet);
    bundle.get(cursor, createSeed);
    bundle.get(cursor, readPhaseFilename);

    // alloy
    bundle.get(cursor, alloyCreateSeed);
    bundle.get(cursor, atom_type::num_atom_types, alloyRatio);

    bundle.get(cursor, collisionStep);
    bundle.get(cursor, 4, collisionLat);
    bundle.get(cursor, pkaEnergy);
    bundle.get(cursor, DIMENSION, direction);

    bundle.get(cursor, potentialFileType);
    bundle.get(cursor, potentialFilename);

    // output section.
    bundle.get(cursor, atomsDumpMode);
    bundle.get(cursor, atomsDumpFilePath);
    bundle.get(cursor, atomsDumpInterval);
    bundle.get(cursor, outByFrame);

    // logs subsection in output section.
    bundle.get(cursor, logs_mode);
    bundle.get(cursor, logs_filename);
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

    os << "simulation.collision.collision_step:" << cv.collisionStep << std::endl;
    os << "simulation.collision.lat:" << cv.collisionLat[0] << "," << cv.collisionLat[1] << ","
       << cv.collisionLat[2] << "," << cv.collisionLat[3] << std::endl;
    os << "pka:" << cv.pkaEnergy << "," ;
    os << "simulation.collision.direction:" << cv.direction[0] << "," << cv.direction[1] << ","
       << cv.direction[2] << std::endl;

    os << "simulation.potential_file.type:" << cv.potentialFileType << std::endl;
    os << "simulation.potential_file.filename:" << cv.potentialFilename << std::endl;

    // output section
    os << "output.mode(copy:0,direct:1):" << cv.atomsDumpMode << std::endl;
    os << "output.dump_interval" << cv.atomsDumpInterval << std::endl;
    os << "output.dump_filename:" << cv.atomsDumpFilePath << std::endl;
    os << "output.logs.mode: " << (cv.logs_mode == LOGS_MODE_CONSOLE ? LOGS_MODE_CONSOLE_STRING : LOGS_MODE_FILE_STRING)
       << "output.logs.by-frame:" << cv.outByFrame << std::endl;
    os << "output.dump_filename:" << cv.logs_filename << std::endl;
    os << "============================================" << std::endl << std::endl;
    return os;
}
