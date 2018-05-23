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
        collisionStep(0), collisionLat{0, 0, 0, 0}, collisionV{0.0, 0.0, 0.0},
        // todo potential type and filename initialize.
        outputMode(OUTPUT_COPY_MODE), outputDumpFilename(DEFAULT_OUTPUT_DUMP_FILENAME) {}

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
    bundle.put(DIMENSION, collisionV);

    bundle.put(potentialFileType);
    bundle.put(potentialFilename);

    // out section
    bundle.put(outputMode);
    bundle.put(outputDumpFilename);
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
    bundle.get(cursor, DIMENSION, collisionV);

    bundle.get(cursor, potentialFileType);
    bundle.get(cursor, potentialFilename);

    bundle.get(cursor, outputMode);
    bundle.get(cursor, outputDumpFilename);
//    }
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
    os << "simulation.collision.collision_v:" << cv.collisionV[0] << "," << cv.collisionV[1] << ","
       << cv.collisionV[2] << std::endl;

    os << "simulation.potential_file.type:" << cv.potentialFileType << std::endl;
    os << "simulation.potential_file.filename:" << cv.potentialFilename << std::endl;

    // output section
    os << "output.mode(copy:0,direct:1):" << cv.outputMode << std::endl;
    os << "output.dump_filename:" << cv.outputDumpFilename << std::endl;
    os << "============================================" << std::endl << std::endl;
    return os;
}