//
// Created by genshen on 2019/10/4.
//

#include <gtest/gtest.h>
#include <atom.h>
#include <world_builder.h>
#include "domain_test_utils.h"

// fixme in mpi, may crash. index out of range
TEST(atom_set_v_test, atom_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    int rand_seek = 1024;
    comm::BccDomain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);
    auto *_atom = new atom(_domain);

    int ra[3] = {97, 0, 0};
    WorldBuilder mWorldBuilder;
    mWorldBuilder.setDomain(_domain)
            .setAtomsContainer(_atom)
            .setBoxSize(space[0], space[1], space[2])
            .setRandomSeed(rand_seek)
            .setLatticeConst(lattice_const)
            .setAlloyRatio(ra)
            .build();

    // test case 1
    const _type_lattice_coord pos[3]{25, 25, 25};
    const double dir[3]{1, 3, 5};
    AtomElement my_atom_before = _atom->atom_list->getAtomEleBySubBoxIndex(2 * pos[0], pos[1], pos[2]);
    _atom->setv(pos, dir, 5000);
    AtomElement my_atom_post = _atom->atom_list->getAtomEleBySubBoxIndex(2 * pos[0], pos[1], pos[2]);

    EXPECT_FLOAT_EQ(my_atom_post.v[0], my_atom_before.v[0] + 222.17975);
    EXPECT_FLOAT_EQ(my_atom_post.v[1], my_atom_before.v[1] + 666.53928);
    EXPECT_FLOAT_EQ(my_atom_post.v[2], my_atom_before.v[2] + 1110.89879);

    // reset
    AtomElement &temp = _atom->atom_list->getAtomEleBySubBoxIndex(2 * pos[0], pos[1], pos[2]);
    temp.v[0] = 0.0;
    temp.v[1] = 0.0;
    temp.v[2] = 0.0;
    // test case 2
    const double dir2[3]{1, 2, 2};
    AtomElement my_atom_before_2 = _atom->atom_list->getAtomEleBySubBoxIndex(2 * pos[0], pos[1], pos[2]);
    _atom->setv(pos, dir2, 20000);
    AtomElement my_atom_post_2 = _atom->atom_list->getAtomEleBySubBoxIndex(2 * pos[0], pos[1], pos[2]);

    EXPECT_FLOAT_EQ(my_atom_post_2.v[0], my_atom_before_2.v[0] + 876.28882);
    EXPECT_FLOAT_EQ(my_atom_post_2.v[1], my_atom_before_2.v[1] + 1752.5776);
    EXPECT_FLOAT_EQ(my_atom_post_2.v[2], my_atom_before_2.v[2] + 1752.5776);
}
