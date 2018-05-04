//
// Created by genshen on 5/2/18.
//

#include <catch2.hpp>
#include <utils/mpi_utils.h>
#include <iostream>
#include "atom.h"
#include "domain.h"

// global domain setting.
int64_t space[3] = {50, 50, 50};
double lattice_const = 0.86;
double cutoff_radius = 1.1421;
int rand_seek = 1024;

DomainDecomposition *getDomainInstance() {
    static DomainDecomposition *_domain = nullptr;
    if (!_domain) {
        return _domain = (new DomainDecomposition(space, lattice_const, cutoff_radius))
                ->decomposition()
                ->createGlobalDomain()
                ->createLocalBoxDomain();
    }
    return _domain;
}

TEST_CASE("domain-test-decomposition", "domain-test") {
    DomainDecomposition *_domain = getDomainInstance();
    atom *_atom = new atom(_domain, lattice_const, cutoff_radius, rand_seek);

    int grid_sum[3] = {0, 0, 0};
    int local_grid[3] = {_domain->getSubBoxLatticeSize(0), _domain->getSubBoxLatticeSize(1),
                         _domain->getSubBoxLatticeSize(2)};
    MPI_Reduce(local_grid, grid_sum, 3, MPI_INT, MPI_SUM, MASTER_PROCESSOR, kiwi::mpiUtils::global_comm);

    if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
        REQUIRE(grid_sum[1] == (2 * space[1])); // todo test
    }
}

void runDomainLocalLatticeSizeTest(int64_t space[3], double lattice_const, double cutoff_radius) {
    DomainDecomposition *_domain = getDomainInstance();

    int nlocalx = floor(_domain->getMeasuredSubBoxUpperBounding(0) / (lattice_const)) -
                  floor(_domain->getMeasuredSubBoxLowerBounding(0) / (lattice_const));
    nlocalx *= 2;
    int nlocaly = floor(_domain->getMeasuredSubBoxUpperBounding(1) / lattice_const) -
                  floor(_domain->getMeasuredSubBoxLowerBounding(1) / lattice_const);
    int nlocalz = floor(_domain->getMeasuredSubBoxUpperBounding(2) / lattice_const) -
                  floor(_domain->getMeasuredSubBoxLowerBounding(2) / lattice_const);

    REQUIRE(nlocalx == _domain->getSubBoxLatticeSize(0));
    REQUIRE(nlocaly == _domain->getSubBoxLatticeSize(1));
    REQUIRE(nlocalz == _domain->getSubBoxLatticeSize(2));
}

TEST_CASE("domain-local-lattice-size", "domain-test") {
    int64_t a[3] = {60, 58, 67};
    runDomainLocalLatticeSizeTest(a, 0.86, 1.1421);
}