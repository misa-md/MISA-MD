//
// Created by genshen(genshenchu@gmail.com) on 2018-3-13.
//

#include "utils/mpi_utils.h"
#include "catch.hpp"
#include "utils/io_writer.h"

TEST_CASE("Write data to shared file", "[MPI_IO_1]") {
//    REQUIRE(Factorial(1) == 1);
    int len = 1024 * 10;
    double *data = new double[len];
    for (int i = 0; i < len; i++) {
        data[i] = mpiUtils::ownRank + 1;
    }
    IOWriter *atomDump = new IOWriter("testfile");
    atomDump->write(data, len);
}

//
//TEST_CASE("Write data to shared file", "[MPI_IO_2]") {
////    REQUIRE(Factorial(1) == 1);
//    char data[1024 * 1123];
//    for (char &i : data) {
//        i = mpiUtils::ownRank;
//    }
//    IOWriter *atomDump = new IOWriter("testfile");
//    atomDump->write(data, 1023 * 1123);
//}

//TEST_CASE("READ_BIN", "BIN") {
//
//}
