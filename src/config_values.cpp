//
// Created by genshen(genshenchu@gmail.com) on 2018-3-9.
//

#include <mpi.h>
#include <iostream>
#include "config_values.h"

ConfigValues::ConfigValues(unsigned int cap) :
        DataPack(cap),
        phaseSpace{0, 0, 0}, cutoffRadius(0.0), latticeConst(0.0),
        timeSteps(10), createPhaseMode(true), createTSet(0.0), createSeed(1),
        readPhaseFilename(""), collisionSteps(0),
        collisionLat{0, 0, 0, 0}, collisionV{0.0, 0.0, 0.0},
        outputMode(OUTPUT_COPY_MODE), outputDumpFilename(DEFAULT_OUTPUT_DUMP_FILENAME) {}

void ConfigValues::packdata() {
    // append data into buffer.
    append(MPI_COMM_WORLD, DIMENSION, phaseSpace); // todo remove MPI_COMM_WORLD to initial method.
    append(MPI_COMM_WORLD, cutoffRadius);
    append(MPI_COMM_WORLD, latticeConst);
    append(MPI_COMM_WORLD, timeSteps);

    append(MPI_COMM_WORLD, createPhaseMode);
    append(MPI_COMM_WORLD, createTSet);
    append(MPI_COMM_WORLD, createSeed);
    append(MPI_COMM_WORLD, readPhaseFilename);

    append(MPI_COMM_WORLD, collisionSteps);
    append(MPI_COMM_WORLD, 4, collisionLat);
    append(MPI_COMM_WORLD, DIMENSION, collisionV);

    append(MPI_COMM_WORLD, potentialFileType);
    append(MPI_COMM_WORLD, potentialFilename);

    // out section
    append(MPI_COMM_WORLD, outputMode);
    append(MPI_COMM_WORLD, outputDumpFilename);
}

void ConfigValues::unpackdata() {
    // fetch data from buffer.
//    if (getPackedData() != nullptr) { // buffer != null
    int cursor = 0;
    recover(MPI_COMM_WORLD, cursor, DIMENSION, phaseSpace);
    recover(MPI_COMM_WORLD, cursor, cutoffRadius);
    recover(MPI_COMM_WORLD, cursor, latticeConst);
    recover(MPI_COMM_WORLD, cursor, timeSteps);

    recover(MPI_COMM_WORLD, cursor, createPhaseMode);
    recover(MPI_COMM_WORLD, cursor, createTSet);
    recover(MPI_COMM_WORLD, cursor, createSeed);
    recover(MPI_COMM_WORLD, cursor, readPhaseFilename);

    recover(MPI_COMM_WORLD, cursor, collisionSteps);
    recover(MPI_COMM_WORLD, cursor, 4, collisionLat);
    recover(MPI_COMM_WORLD, cursor, DIMENSION, collisionV);

    recover(MPI_COMM_WORLD, cursor, potentialFileType);
    recover(MPI_COMM_WORLD, cursor, potentialFilename);

    recover(MPI_COMM_WORLD, cursor, outputMode);
    recover(MPI_COMM_WORLD, cursor, outputDumpFilename);
//    }
}

std::ostream &operator<<(std::ostream &os, const ConfigValues &cv) {
    // simulation section
    os << "===========config of simulation=============" << std::endl;
    os << "simulation.phase_space:" << cv.phaseSpace[0] << "," << cv.phaseSpace[1]
       << "," << cv.phaseSpace[2] << "," << endl;
    os << "simulation.cutoff_radius:" << cv.cutoffRadius << endl;
    os << "simulation.lattice_const:" << cv.latticeConst << endl;
    os << "simulation.timesteps:" << cv.timeSteps << endl;

    os << "simulation.createphase.createPhaseMode:" << (cv.createPhaseMode ? "true" : "false") << endl;
    os << "simulation.createphase.tSet:" << cv.createTSet << endl;
    os << "simulation.createphase.seed:" << cv.createSeed << endl;
    os << "simulation.createphase.readPhaseFilename:" << cv.readPhaseFilename << endl;

    os << "simulation.collision.collision_steps:" << cv.collisionSteps << endl;
    os << "simulation.collision.lat:" << cv.collisionLat[0] << "," << cv.collisionLat[1] << ","
       << cv.collisionLat[2] << "," << cv.collisionLat[3] << endl;
    os << "simulation.collision.collision_v:" << cv.collisionV[0] << "," << cv.collisionV[1] << ","
       << cv.collisionV[2] << endl;

    os << "simulation.potential_file.type:" << cv.potentialFileType << endl;
    os << "simulation.potential_file.filename:" << cv.potentialFilename << endl;

    // output section
    os << "output.mode(copy:0,direct:1):" << cv.outputMode << endl;
    os << "output.dump_filename:" << cv.outputDumpFilename << endl;
    os << "============================================" << std::endl << std::endl;
    return os;
}
