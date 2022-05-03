//
// Created by genshen on 2018-12-31.
//

#include <cmath>
#include <types/pre_define.h>
#include <comm/domain/domain.h>

#include "atom_set.h"
#include "lattice/ws_utils.h"

AtomSet::AtomSet(const double cutoff_radius,
                 const _type_lattice_size extended_lattice_size[DIMENSION],
                 const _type_lattice_size sub_box_lattice_size[DIMENSION],
                 const _type_lattice_size ghost_lattice_size[DIMENSION])
        : _cutoffRadius(cutoff_radius) {
    // the length of atom array at x direction is doubled due to the special data structure.
    atom_list = new AtomList(extended_lattice_size[0] * 2,
                             extended_lattice_size[1],
                             extended_lattice_size[2],
                             sub_box_lattice_size[0] * 2,
                             sub_box_lattice_size[1],
                             sub_box_lattice_size[2],
                             ghost_lattice_size[0] * 2,
                             ghost_lattice_size[1],
                             ghost_lattice_size[2]);

    inter_atom_list = new InterAtomList();
    // create neighbour relative index.
    neighbours = new NeighbourIndex<_type_neighbour_index_ele>(atom_list->_atoms._data(), atom_list->lattice);
}

AtomSet::~AtomSet() {
    delete atom_list;
    delete inter_atom_list;
    delete neighbours;
}

void AtomSet::calcNeighbourIndices(const double cutoff_radius_factor, const _type_lattice_size cut_lattice) {
    neighbours->make(cut_lattice, cutoff_radius_factor);
}

void
AtomSet::addAtom(comm::BccDomain *p_domain, const AtomElement atom) {
    if (ws::isInBox(atom.x[0], atom.x[1], atom.x[2], p_domain)) {
        const _type_atom_index near_atom_index = ws::findNearLatIndexInSubBox(atom_list->lattice, atom, p_domain);
        MD_LOAD_ATOM_VAR(atom_near, atom_list, near_atom_index);
        if (MD_GET_ATOM_TYPE(atom_near, near_atom_index) == atom_type::INVALID) {
            MD_SET_ATOM_ID(atom_near, near_atom_index, atom.id);
            MD_SET_ATOM_TYPE(atom_near, near_atom_index, atom.type);

            MD_SET_ATOM_X(atom_near, near_atom_index, 0, atom.x[0]);
            MD_SET_ATOM_X(atom_near, near_atom_index, 1, atom.x[1]);
            MD_SET_ATOM_X(atom_near, near_atom_index, 2, atom.x[2]);

            MD_SET_ATOM_V(atom_near, near_atom_index, 0, atom.v[0]);
            MD_SET_ATOM_V(atom_near, near_atom_index, 1, atom.v[1]);
            MD_SET_ATOM_V(atom_near, near_atom_index, 2, atom.v[2]);
        }else {
            inter_atom_list->addInterAtom(atom);
        }
    }
}

_type_atom_count AtomSet::getnlocalatom(comm::Domain *p_domain) {
    return (p_domain->sub_box_lattice_size[0] * p_domain->sub_box_lattice_size[1] *
            p_domain->sub_box_lattice_size[2]);
}

#ifdef MD_RUNTIME_CHECKING

_type_atom_count AtomSet::realAtoms() const {
    _type_atom_count count = 0;
    atom_list->foreachSubBoxAtom(
            [=, &count](const _type_atom_index gid) {
                MD_LOAD_ATOM_VAR(_atom_ref, atom_list, gid);
                if (!MD_IS_ATOM_TYPE_INTER(_atom_ref, gid)) {
                    count++;
                }
            }
    );
    return count;
}

#endif
