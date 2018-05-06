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

Domain *getDomainInstance() {
    static Domain *_domain = nullptr;
    if (!_domain) {
        return _domain = (new Domain(space, lattice_const, cutoff_radius))
                ->decomposition()
                ->createGlobalDomain()
                ->createSubBoxDomain();
    }
    return _domain;
}

TEST_CASE("domain-test-decomposition", "domain-test") {
    Domain *_domain = getDomainInstance();
    atom *_atom = new atom(_domain, lattice_const, cutoff_radius, rand_seek);

    int grid_sum[3] = {0, 0, 0};
    int local_grid[3] = {_domain->getSubBoxLatticeSize(0), _domain->getSubBoxLatticeSize(1),
                         _domain->getSubBoxLatticeSize(2)};
    MPI_Reduce(local_grid, grid_sum, 3, MPI_INT, MPI_SUM, MASTER_PROCESSOR, kiwi::mpiUtils::global_comm);

    if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
        REQUIRE(grid_sum[1] == (2 * space[1])); // todo test
    }
}

TEST_CASE("domain-local-lattice-size", "domain-test") {
    Domain *_domain = getDomainInstance();

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

TEST_CASE("domain-ghost-lattice-size", "domain-test") {
    Domain *_domain = getDomainInstance();

    int nghostx = _domain->getSubBoxLatticeSize(0) + 2 * 2 * ceil(_domain->getMeasuredGhostLength(0) / lattice_const);
    int nghosty = _domain->getSubBoxLatticeSize(1) + 2 * ceil(_domain->getMeasuredGhostLength(1) / lattice_const);
    int nghostz = _domain->getSubBoxLatticeSize(2) + 2 * ceil(_domain->getMeasuredGhostLength(2) / lattice_const);

    REQUIRE(nghostx == _domain->getGhostLatticeSize(0));
    REQUIRE(nghosty == _domain->getGhostLatticeSize(1));
    REQUIRE(nghostz == _domain->getGhostLatticeSize(2));
}

TEST_CASE("domain-local-lattice-coord", "domain-test") {
    Domain *_domain = getDomainInstance();

    // lower boundary of lattice coordinate of local sub-box
    int lolocalx = floor(_domain->getMeasuredSubBoxLowerBounding(0) / lattice_const) * 2;
    int lolocaly = floor(_domain->getMeasuredSubBoxLowerBounding(1) / lattice_const);
    int lolocalz = floor(_domain->getMeasuredSubBoxLowerBounding(2) / lattice_const);

    REQUIRE(lolocalx == _domain->getSubBoxLatticeCoordLower(0));
    REQUIRE(lolocaly == _domain->getSubBoxLatticeCoordLower(1));
    REQUIRE(lolocalz == _domain->getSubBoxLatticeCoordLower(2));

    // upper boundary of lattice coordinate of local sub-box
    int uplocalx = floor(_domain->getMeasuredSubBoxUpperBounding(0) / lattice_const) * 2;
    int uplocaly = floor(_domain->getMeasuredSubBoxUpperBounding(1) / lattice_const);
    int uplocalz = floor(_domain->getMeasuredSubBoxUpperBounding(2) / lattice_const);

    REQUIRE(uplocalx == _domain->getSubBoxLatticeCoordUpper(0));
    REQUIRE(uplocaly == _domain->getSubBoxLatticeCoordUpper(1));
    REQUIRE(uplocalz == _domain->getSubBoxLatticeCoordUpper(2));
}

TEST_CASE("domain-ghost-lattice-coord", "domain-test") {
    Domain *_domain = getDomainInstance();

    // lower boundary of lattice coordinate of ghost
    int loghostx = _domain->getSubBoxLatticeCoordLower(0) - 2 * ceil(cutoff_radius / lattice_const);
    int loghosty = _domain->getSubBoxLatticeCoordLower(1) - ceil(cutoff_radius / lattice_const);
    int loghostz = _domain->getSubBoxLatticeCoordLower(2) - ceil(cutoff_radius / lattice_const);

    REQUIRE(loghostx == _domain->getGhostLatticeCoordLower(0));
    REQUIRE(loghosty == _domain->getGhostLatticeCoordLower(1));
    REQUIRE(loghostz == _domain->getGhostLatticeCoordLower(2));

    // upper boundary of lattice coordinate of ghost
    int upghostx = _domain->getSubBoxLatticeCoordUpper(0) + 2 * ceil(cutoff_radius / lattice_const);
    int upghosty = _domain->getSubBoxLatticeCoordUpper(1) + ceil(cutoff_radius / lattice_const);
    int upghostz = _domain->getSubBoxLatticeCoordUpper(2) + ceil(cutoff_radius / lattice_const);

    REQUIRE(upghostx == _domain->getGhostLatticeCoordUpper(0));
    REQUIRE(upghosty == _domain->getGhostLatticeCoordUpper(1));
    REQUIRE(upghostz == _domain->getGhostLatticeCoordUpper(2));
}
