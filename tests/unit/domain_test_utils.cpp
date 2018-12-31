//
// Created by genshen on 5/12/18.
//

#include <utils/mpi_utils.h>
#include "domain_test_utils.h"

Domain *getDomainInstance(int64_t space[3], double lattice_const, double cutoff_radius) {
    MPI_Comm mpi_comm;
    Domain *p_domain = Domain::Builder()
            .setComm(kiwi::mpiUtils::global_process, &mpi_comm)
            .setPhaseSpace(space)
            .setCutoffRadius(cutoff_radius)
            .setLatticeConst(lattice_const)
            .build();
    kiwi::mpiUtils::onGlobalCommChanged(mpi_comm); // notice MPI changed
    return p_domain;
}
