//
// Created by genshen on 2019/10/3.
//

#include <gtest/gtest.h>
#include <system_configuration.h>

#include "fixtures/world_builder_test_fixture.h"

TEST_F(WorldBuilderTestFixture, configuration_temperature_test) {
    configuration::rescale(600, static_cast<_type_atom_count>(2 * space[0] * space[1] * space[2]),
                           _atom->atom_list, _atom->inter_atom_list);
    // test temperature
    const double T = configuration::temperature(static_cast<_type_atom_count>(2 * space[0] * space[1] * space[2]),
                                                _atom->atom_list, _atom->inter_atom_list);
    EXPECT_FLOAT_EQ(T, 600);

    // test configuration system temperature.
    const double e = configuration::kineticEnergy(_atom->getAtomList(), _atom->getInterList(),
                                                  configuration::ReturnMod::All, 0);
    const double T2 = configuration::temperature(e, 2 * space[0] * space[1] * space[2]);
    EXPECT_DOUBLE_EQ(T, T2);
}
