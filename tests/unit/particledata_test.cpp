//
// Created by genshen on 2018-12-02.
//

#include <gtest/gtest.h>
#include <pack/particledata.h>
#include <utils/mpi_utils.h>

TEST(particledata_test, particledata_test_asny) {
    MPI_Datatype _mpi_particle_data;
    particledata::setMPIType(_mpi_particle_data);

    kiwi::RID rank = kiwi::mpiUtils::global_process.own_rank;
    kiwi::RID allrank = kiwi::mpiUtils::global_process.all_ranks;

    particledata p;
    p.id = 12;
    p.r[0] = 8.0;

    // send
    MPI_Request send_requests;
    MPI_Isend(&p, 1, _mpi_particle_data, (rank + 1) % allrank, 0x100, MPI_COMM_WORLD, &send_requests);
    // receive
    particledata p2;
    MPI_Request rec_requests;
    MPI_Irecv(&p2, 1, _mpi_particle_data, (rank + allrank - 1) % allrank, 0x100,
              MPI_COMM_WORLD, &rec_requests);
    // wait
    MPI_Status send_status, rec_status;
    MPI_Wait(&send_requests, &send_status);
    MPI_Wait(&rec_requests, &rec_status);

    EXPECT_EQ(p2.id, 12);
    EXPECT_DOUBLE_EQ(p2.r[0], 8.0); // fixme  Which is: 0, To be equal to: 8.0
}

TEST(particledata_test, particledata_test_sync_reci) {
    MPI_Datatype _mpi_particle_data;
    particledata::setMPIType(_mpi_particle_data);

    kiwi::RID rank = kiwi::mpiUtils::global_process.own_rank;
    kiwi::RID allrank = kiwi::mpiUtils::global_process.all_ranks;

    particledata p;
    p.id = 12;
    p.r[0] = 8.0;

    MPI_Status send_status, rec_status;
    // send
    MPI_Request send_requests;
    MPI_Isend(&p, 1, _mpi_particle_data, (rank + 1) % allrank, 0x100, MPI_COMM_WORLD, &send_requests);
    // receive
    particledata p2;
    MPI_Recv(&p2, 1, _mpi_particle_data, (rank + allrank - 1) % allrank, 0x100,
             MPI_COMM_WORLD, &rec_status);
    // wait
    MPI_Wait(&send_requests, &send_status);

    EXPECT_EQ(p2.id, 12);
    EXPECT_DOUBLE_EQ(p2.r[0], 8.0); // fixme  Which is: 0, To be equal to: 8.0
}
