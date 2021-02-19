//
// Created by genshen on 2021/2/19.
//

#include "world_builder_test_fixture.h"

void WorldBuilderTestFixture::SetUp() {
    DomainFixture::SetUp();

    _atom = new atom(p_domain);
    mWorldBuilder.setDomain(p_domain)
            .setAtomsContainer(_atom)
            .setBoxSize(space[0], space[1], space[2])
            .setRandomSeed(rand_seek)
            .setLatticeConst(lattice_const)
            .setAlloyRatio(ra)
            .build();
}
