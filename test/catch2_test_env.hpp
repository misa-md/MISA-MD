//
// Created by genshen on 2018-3-12.
//

#ifndef CRYSTAL_MD_CATCH2_TEST_ENV_H
#define CRYSTAL_MD_CATCH2_TEST_ENV_H


#define CATCH_CONFIG_COLOUR_ANSI
#define CATCH_CONFIG_RUNNER
// #define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file

#include <iostream>
#include "test_config.h"
#include "catch.hpp"

#ifdef TEST_MPI_ENABLE_FLAG

#include "mpi.h"

#endif  // end TEST_MPI_ENABLE_FLAG


#endif //CRYSTAL_MD_CATCH2_TEST_ENV_H
