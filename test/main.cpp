//
// Created by genshen on 2018-3-12.
//

#include "catch2_test_env.hpp"

int main(int argc, char *argv[]) {
    Catch::Session session;
    // global setup...
#ifdef TEST_MPI_ENABLE_FLAG
    // Initialize the MPI environment.
    MPI_Init(&argc, &argv);
    std::cout << "initialed mpi." << std::endl;
#endif  // end TEST_MPI_ENABLE_FLAG

    int result = session.run(argc, argv);

// global clean-up...
#ifdef TEST_MPI_ENABLE_FLAG
    // Finalize the MPI environment.
    MPI_Finalize();
#endif  // end TEST_MPI_ENABLE_FLAG

    return result;
}
