//
// Created by genshen on 2019/10/3.
//

#include <gtest/gtest.h>
#include <system_configuration.h>
#include <atom.h>
#include <world_builder.h>

#include "domain_test_utils.h"

TEST(configuration_temperature_test, configuration_test) {
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
