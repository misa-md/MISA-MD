//
// Created by genshen on 2018-05-02.
//

#include <gtest/gtest.h>
#include <utils/mpi_utils.h>
#include <iostream>
#include <cmath>
#include "domain_test_utils.h"
#include "atom.h"

TEST(domain_test_decomposition, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    int rand_seek = 1024;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);
    auto *_atom = new atom(_domain, lattice_const, cutoff_radius_factor, rand_seek);

    int grid_sum[3] = {0, 0, 0};
    int local_grid[3] = {_domain->getSubBoxLatticeSize(0), _domain->getSubBoxLatticeSize(1),
                         _domain->getSubBoxLatticeSize(2)};
    MPI_Reduce(local_grid, grid_sum, 3, MPI_INT, MPI_SUM, MASTER_PROCESSOR, kiwi::mpiUtils::global_comm);

    if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
        EXPECT_EQ(grid_sum[1], (2 * space[1])); // todo test
    }

    delete _domain;
    delete _atom;
}

TEST(domain_local_lattice_size, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    int nlocalx = floor(_domain->getMeasuredSubBoxUpperBounding(0) / (lattice_const)) -
                  floor(_domain->getMeasuredSubBoxLowerBounding(0) / (lattice_const));
    nlocalx *= 2;
    int nlocaly = floor(_domain->getMeasuredSubBoxUpperBounding(1) / lattice_const) -
                  floor(_domain->getMeasuredSubBoxLowerBounding(1) / lattice_const);
    int nlocalz = floor(_domain->getMeasuredSubBoxUpperBounding(2) / lattice_const) -
                  floor(_domain->getMeasuredSubBoxLowerBounding(2) / lattice_const);

    EXPECT_EQ(nlocalx, _domain->getSubBoxLatticeSize(0));
    EXPECT_EQ(nlocaly, _domain->getSubBoxLatticeSize(1));
    EXPECT_EQ(nlocalz, _domain->getSubBoxLatticeSize(2));

    delete _domain;
}

TEST(domain_ghost_lattice_size, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    int nghostx = _domain->getSubBoxLatticeSize(0) + 2 * 2 * ceil(_domain->getMeasuredGhostLength(0) / lattice_const);
    int nghosty = _domain->getSubBoxLatticeSize(1) + 2 * ceil(_domain->getMeasuredGhostLength(1) / lattice_const);
    int nghostz = _domain->getSubBoxLatticeSize(2) + 2 * ceil(_domain->getMeasuredGhostLength(2) / lattice_const);

    EXPECT_EQ(nghostx, _domain->getGhostExtLatticeSize(0));
    EXPECT_EQ(nghosty, _domain->getGhostExtLatticeSize(1));
    EXPECT_EQ(nghostz, _domain->getGhostExtLatticeSize(2));

    delete _domain;
}

TEST(domain_local_lattice_coord, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    // lower boundary of lattice coordinate of local sub-box
    int lolocalx = floor(_domain->getMeasuredSubBoxLowerBounding(0) / lattice_const) * 2;
    int lolocaly = floor(_domain->getMeasuredSubBoxLowerBounding(1) / lattice_const);
    int lolocalz = floor(_domain->getMeasuredSubBoxLowerBounding(2) / lattice_const);

    EXPECT_EQ(lolocalx, _domain->getGlobalSubBoxLatticeCoordLower(0));
    EXPECT_EQ(lolocaly, _domain->getGlobalSubBoxLatticeCoordLower(1));
    EXPECT_EQ(lolocalz, _domain->getGlobalSubBoxLatticeCoordLower(2));

    // upper boundary of lattice coordinate of local sub-box
    int uplocalx = floor(_domain->getMeasuredSubBoxUpperBounding(0) / lattice_const) * 2;
    int uplocaly = floor(_domain->getMeasuredSubBoxUpperBounding(1) / lattice_const);
    int uplocalz = floor(_domain->getMeasuredSubBoxUpperBounding(2) / lattice_const);

    EXPECT_EQ(uplocalx, _domain->getGlobalSubBoxLatticeCoordUpper(0));
    EXPECT_EQ(uplocaly, _domain->getGlobalSubBoxLatticeCoordUpper(1));
    EXPECT_EQ(uplocalz, _domain->getGlobalSubBoxLatticeCoordUpper(2));

    delete _domain;
}

TEST(domain_ghost_lattice_coord, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    // lower boundary of lattice coordinate of ghost
    int loghostx = _domain->getGlobalSubBoxLatticeCoordLower(0) - 2 * ceil(cutoff_radius_factor / lattice_const);
    int loghosty = _domain->getGlobalSubBoxLatticeCoordLower(1) - ceil(cutoff_radius_factor / lattice_const);
    int loghostz = _domain->getGlobalSubBoxLatticeCoordLower(2) - ceil(cutoff_radius_factor / lattice_const);

    EXPECT_EQ(loghostx, _domain->getGlobalGhostLatticeCoordLower(0));
    EXPECT_EQ(loghosty, _domain->getGlobalGhostLatticeCoordLower(1));
    EXPECT_EQ(loghostz, _domain->getGlobalGhostLatticeCoordLower(2));

    // upper boundary of lattice coordinate of ghost
    int upghostx = _domain->getGlobalSubBoxLatticeCoordUpper(0) + 2 * ceil(cutoff_radius_factor / lattice_const);
    int upghosty = _domain->getGlobalSubBoxLatticeCoordUpper(1) + ceil(cutoff_radius_factor / lattice_const);
    int upghostz = _domain->getGlobalSubBoxLatticeCoordUpper(2) + ceil(cutoff_radius_factor / lattice_const);

    EXPECT_EQ(upghostx, _domain->getGlobalGhostLatticeCoordUpper(0));
    EXPECT_EQ(upghosty, _domain->getGlobalGhostLatticeCoordUpper(1));
    EXPECT_EQ(upghostz, _domain->getGlobalGhostLatticeCoordUpper(2));

    delete _domain;
}
