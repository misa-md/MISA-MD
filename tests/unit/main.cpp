//
// Created by genshen on 2018-3-12.
//

#include <test/gtest_env.h>
#include <utils/mpi_domain.h>
#include "test_config.h"

// see https://github.com/google/googletest/issues/822 for more information.
// main function for adapting mpi environment
int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef TEST_MPI_ENABLE_FLAG
    ::testing::AddGlobalTestEnvironment(new kiwi::MPIEnvironment);
#endif  // end TEST_MPI_ENABLE_FLAG
    return RUN_ALL_TESTS();
}
