//
// Created by genshen on 5/12/18.
//

#include <gtest/gtest.h>
#include <world_builder.h>
#include <iostream>

#include "domain_test_utils.h"

TEST(zero_momentum_test, world_builder_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    int rand_seek = 1024;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);
    auto *_atom = new atom(_domain, lattice_const, cutoff_radius_factor, rand_seek);

    int ra[3] = {97, 3, 1};
    WorldBuilder mWorldBuilder;
    mWorldBuilder.setDomain(_domain)
            .setAtomsContainer(_atom)
            .setBoxSize(space[0], space[1], space[2])
            .setRandomSeed(rand_seek)
            .setLatticeConst(lattice_const)
            .setTset(600)
            .setAlloyRatio(ra)
            .build();

    double p[4];
    mWorldBuilder.vcm(p);
    std::cout << "<<<<<"
              << mWorldBuilder.computeScalar(static_cast<_type_atom_count>(2 * space[0] * space[1] * space[2]));

    for (int i = 0; i < DIMENSION; i++) {
        double a = p[i];
//        double b = a + std::numeric_limits<double>::epsilon();
//        REQUIRE_FALSE(b == a);
        EXPECT_FLOAT_EQ(0, a);  // zero momentumat each dimension.
    }
}
