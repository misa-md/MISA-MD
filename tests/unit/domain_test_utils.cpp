//
// Created by genshen on 5/12/18.
//

#include <utils/mpi_utils.h>
#include "domain_test_utils.h"

comm::Domain *getDomainInstance(int64_t space[3], double lattice_const, double cutoff_radius) {
    MPI_Comm mpi_comm; // new domain
    comm::mpi_process m_process{};
    m_process.comm = MPI_COMM_WORLD;
    MPI_Comm_size(m_process.comm, &m_process.all_ranks);
    MPI_Comm_rank(m_process.comm, &m_process.own_rank);

    comm::Domain *p_domain = comm::Domain::Builder()
            .setComm(m_process, &mpi_comm)
            .setPhaseSpace(space)
            .setCutoffRadius(cutoff_radius)
            .setLatticeConst(lattice_const)
            .build();
    kiwi::mpiUtils::onGlobalCommChanged(mpi_comm); // notice MPI changed
    return p_domain;
}
