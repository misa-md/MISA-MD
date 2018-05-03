//
// Created by genshen on 5/2/18.
//

#include <catch2.hpp>
#include <utils/mpi_utils.h>
#include <iostream>
#include "atom.h"
#include "domain.h"

TEST_CASE("domain-test-decomposition", "domain-test") {
    int64_t space[3] = {50, 50, 50};
    double lattice_const = 0.86;
    double cutoff_radius = 1.1421;
    int rand_seek = 1024;

    DomainDecomposition *_domain_decomposition = (new DomainDecomposition())
            ->decomposition()
            ->createGlobalDomain(space, lattice_const)
            ->createLocalBoxDomain(space, lattice_const, cutoff_radius);
    atom *_atom = new atom(_domain_decomposition, lattice_const, cutoff_radius, rand_seek);

    int grid_sum[3] = {0, 0, 0};
    int local_grid[3] = {_atom->nlocalx, _atom->nlocaly, _atom->nlocalz};
    MPI_Reduce(local_grid, grid_sum, 3, MPI_INT, MPI_SUM, MASTER_PROCESSOR, kiwi::mpiUtils::global_comm);

    if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
        REQUIRE(grid_sum[1] == (2 * space[1])); // todo test
//        std::cout << "local_grid " << local_grid[0] << "|" << local_grid[1] << "|" << local_grid[2];
    }
}
