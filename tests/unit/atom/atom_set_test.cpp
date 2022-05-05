//
// Created by genshen on 2022/5/3.
//

#include <gtest/gtest.h>
#include <atom/neighbour_index.h>

#include "../fixtures/domain_test_fixture.h"
#include "atom.h"

TEST_F(DomainFixture, atom_set_add_atom_test) {
    atom *_atom_container = new atom(p_domain);

    // set all types to invalid
    _atom_container->atom_list->foreachSubBoxAtom([_atom_container](_type_atom_index gid) {
        MD_LOAD_ATOM_VAR(_atom_ref, _atom_container->atom_list, gid);
        MD_SET_ATOM_TYPE(_atom_ref, gid, atom_type::INVALID);
    });

    // add atoms into atom container
    AtomElement atom{};
    atom.id = 10;
    atom.x[0] = lattice_const * 10;
    atom.x[1] = lattice_const * 10;
    atom.x[2] = lattice_const * 10;
    _atom_container->addAtom(p_domain, atom);
    // check
    EXPECT_EQ(_atom_container->atom_list->countValidAtoms(), 1);

    // add another atom and check.
    atom.id = 11;
    _atom_container->addAtom(p_domain, atom);
    EXPECT_EQ(_atom_container->atom_list->countValidAtoms(), 1);
    EXPECT_EQ(_atom_container->inter_atom_list->nLocalInter(), 1);
    EXPECT_EQ(_atom_container->inter_atom_list->inter_list.front().id, 11);

    delete _atom_container;
}
