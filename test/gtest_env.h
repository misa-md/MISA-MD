//
// Created by genshen on 6/2/18.
//

#ifndef CRYSTALMD_GTEST_ENV_H
#define CRYSTALMD_GTEST_ENV_H

#include <gtest/gtest.h>
#include <utils/mpi_utils.h>
#include "test_config.h"

#ifdef TEST_MPI_ENABLE_FLAG

#include "mpi.h"

class MPIEnvironment : public ::testing::Environment {
public:
    void SetUp() override {
        char **argv = nullptr;
        int argc = 0;
        kiwi::mpiUtils::initialMPI(argc, argv);
    }

    void TearDown() override {
        kiwi::mpiUtils::finishMPI();
    }

    ~MPIEnvironment() override = default;
};

#endif  // end TEST_MPI_ENABLE_FLAG


#endif //CRYSTALMD_GTEST_ENV_H
