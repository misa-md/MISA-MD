//
// Created by genshen on 6/2/18.
//

#ifndef CRYSTALMD_GTEST_ENV_H
#define CRYSTALMD_GTEST_ENV_H

#include <gtest/gtest.h>
#include "test_config.h"

#ifdef TEST_MPI_ENABLE_FLAG

#include "mpi.h"

class MPIEnvironment : public ::testing::Environment {
public:
    void SetUp() override {
        char **argv;
        int argc = 0;
        int mpiError = MPI_Init(&argc, &argv);
        ASSERT_FALSE(mpiError);
    }

    void TearDown() override {
        int mpiError = MPI_Finalize();
        ASSERT_FALSE(mpiError);
    }

    ~MPIEnvironment() override = default;
};

#endif  // end TEST_MPI_ENABLE_FLAG


#endif //CRYSTALMD_GTEST_ENV_H
