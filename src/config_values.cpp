//
// Created by genshen(genshenchu@gmail.com) on 2018-3-9.
//

#include <mpi.h>
#include <iostream>
#include <utils/bundle.h>
#include "config_values.h"

ConfigValues::ConfigValues() :
        phaseSpace{0, 0, 0}, cutoffRadiusFactor(0.0), latticeConst(0.0),
        timeSteps(10), createPhaseMode(true), createTSet(0.0), createSeed(1),
        readPhaseFilename(""), collisionSteps(0),
        collisionLat{0, 0, 0, 0}, collisionV{0.0, 0.0, 0.0},
        outputMode(OUTPUT_COPY_MODE), outputDumpFilename(DEFAULT_OUTPUT_DUMP_FILENAME) {}

void ConfigValues::packdata(kiwi::Bundle &bundle) {
    // append data into buffer.
    bundle.put(MPI_COMM_WORLD, DIMENSION, phaseSpace); // todo remove MPI_COMM_WORLD to initial method.
    bundle.put(MPI_COMM_WORLD, cutoffRadiusFactor);
    bundle.put(MPI_COMM_WORLD, latticeConst);
    bundle.put(MPI_COMM_WORLD, timeSteps);

    bundle.put(MPI_COMM_WORLD, createPhaseMode);
    bundle.put(MPI_COMM_WORLD, createTSet);
    bundle.put(MPI_COMM_WORLD, createSeed);
    bundle.put(MPI_COMM_WORLD, readPhaseFilename);

    bundle.put(MPI_COMM_WORLD, collisionSteps);
    bundle.put(MPI_COMM_WORLD, 4, collisionLat);
    bundle.put(MPI_COMM_WORLD, DIMENSION, collisionV);

    bundle.put(MPI_COMM_WORLD, potentialFileType);
    bundle.put(MPI_COMM_WORLD, potentialFilename);

    // out section
    bundle.put(MPI_COMM_WORLD, outputMode);
    bundle.put(MPI_COMM_WORLD, outputDumpFilename);
}

void ConfigValues::unpackdata(kiwi::Bundle &bundle) {
    // fetch data from buffer.
//    if (getPackedData() != nullptr) { // buffer != null
    int cursor = 0;
    bundle.get(MPI_COMM_WORLD, cursor, DIMENSION, phaseSpace);
    bundle.get(MPI_COMM_WORLD, cursor, cutoffRadiusFactor);
    bundle.get(MPI_COMM_WORLD, cursor, latticeConst);
    bundle.get(MPI_COMM_WORLD, cursor, timeSteps);

    bundle.get(MPI_COMM_WORLD, cursor, createPhaseMode);
    bundle.get(MPI_COMM_WORLD, cursor, createTSet);
    bundle.get(MPI_COMM_WORLD, cursor, createSeed);
    bundle.get(MPI_COMM_WORLD, cursor, readPhaseFilename);

    bundle.get(MPI_COMM_WORLD, cursor, collisionSteps);
    bundle.get(MPI_COMM_WORLD, cursor, 4, collisionLat);
    bundle.get(MPI_COMM_WORLD, cursor, DIMENSION, collisionV);

    bundle.get(MPI_COMM_WORLD, cursor, potentialFileType);
    bundle.get(MPI_COMM_WORLD, cursor, potentialFilename);

    bundle.get(MPI_COMM_WORLD, cursor, outputMode);
    bundle.get(MPI_COMM_WORLD, cursor, outputDumpFilename);
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

    os << "simulation.collision.collision_steps:" << cv.collisionSteps << std::endl;
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
