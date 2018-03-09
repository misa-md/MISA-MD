//
// Created by genshen(genshenchu@gmail.com) on 2017/5/7.
//
#include "mpi.h"

#ifndef CRYSTAL_MD_MPI_UTILS_H
#define CRYSTAL_MD_MPI_UTILS_H

#define MASTER_PROCESSOR 0

namespace mpiUtils {
    // you can use ownRank and allRanks after called function initialMPI in tiny_fmm.h
    extern int ownRank;
    extern int allRanks;

    void initialMPI();

    void finishMPI();

};

#endif // CRYSTAL_MD_MPI_UTILS_H
