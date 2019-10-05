//
// Created by genshen on 5/12/18.
//

#include <gtest/gtest.h>
#include <world_builder.h>
#include <system_configuration.h>
#include <utils/mpi_utils.h>

#include "domain_test_utils.h"

// @MPI
TEST(zero_momentum_test, world_builder_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    int rand_seek = 1024;
    comm::BccDomain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);
    auto *_atom = new atom(_domain);

    int ra[3] = {97, 3, 1};
    WorldBuilder mWorldBuilder;
    mWorldBuilder.setDomain(_domain)
            .setAtomsContainer(_atom)
            .setBoxSize(space[0], space[1], space[2])
            .setRandomSeed(rand_seek)
            .setLatticeConst(lattice_const)
            .setAlloyRatio(ra)
            .build();

    // test zero momentum.
    double p[4];
    mWorldBuilder.vcm(p);
    for (int i = 0; i < DIMENSION; i++) {
        EXPECT_NEAR(0, p[i], 1e-8);  // zero momentum in each dimension.
    }
}
