//
// Created by genshen on 5/12/18.
//

#include <gtest/gtest.h>
#include <world_builder.h>
#include <system_configuration.h>
#include <utils/mpi_utils.h>
#include "fixtures/domain_test_fixture.h"

// @MPI
TEST_F(DomainFixture, world_build_zero_momentum_test) {
    int rand_seek = 1024;
    comm::BccDomain *_domain = p_domain;
    auto *_atom = new atom(_domain);

    std::vector<int> ra = {97, 3, 1};
    WorldBuilder mWorldBuilder;
    mWorldBuilder.setDomain(_domain)
            .setAtomsContainer(_atom)
            .setBoxSize(space[0], space[1], space[2])
            .setRandomSeed(rand_seek)
            .setLatticeConst(lattice_const)
            .setTset(600)
            .setAlloyRatio(ra)
            .build();

    // test configuration system temperature.
    const double e = configuration::kineticEnergy(_atom->getAtomList(), _atom->getInterList(),
                                                  configuration::ReturnMod::All, 0);
    const double T2 = configuration::temperature(e, 2 * space[0] * space[1] * space[2]);
    EXPECT_FLOAT_EQ(T2, 600);

    // test zero momentum.
    double p[4];
    mWorldBuilder.vcm(p);
    for (int i = 0; i < DIMENSION; i++) {
        EXPECT_NEAR(0, p[i], 1e-8);  // zero momentum in each dimension.
    }
}
