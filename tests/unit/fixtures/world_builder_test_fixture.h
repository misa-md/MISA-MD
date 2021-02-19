//
// Created by genshen on 2021/2/19.
//

#ifndef MISA_MD_WORLD_BUILDER_TEST_FIXTURE_H
#define MISA_MD_WORLD_BUILDER_TEST_FIXTURE_H

#include <atom.h>
#include <world_builder.h>
#include "domain_test_fixture.h"

class WorldBuilderTestFixture : public DomainFixture {
protected:
    void SetUp() override;

protected:
    static constexpr int rand_seek = 1024;
    atom *_atom = nullptr;
    int ra[3] = {97, 3, 1};
    WorldBuilder mWorldBuilder;
};


#endif //MISA_MD_WORLD_BUILDER_TEST_FIXTURE_H
