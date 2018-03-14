//
// Created by genshen on 2018-3-12.
//

#include "catch2_test_env.hpp"

int main(int argc, char *argv[]) {
    Catch::Session session;
    // global setup...
#ifdef TEST_MPI_ENABLE_FLAG
    // Initialize the MPI environment.
    mpiUtils::initialMPI();
#endif  // end TEST_MPI_ENABLE_FLAG

    int result = session.run(argc, argv);

// global clean-up...
#ifdef TEST_MPI_ENABLE_FLAG
    // Finalize the MPI environment.
    mpiUtils::finishMPI();
#endif  // end TEST_MPI_ENABLE_FLAG

    return result;
}
