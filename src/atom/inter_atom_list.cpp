//
// Created by genshen on 5/19/18.
//

#include <cmath>
#include <comm/domain/domain.h>
#include <comm/comm.hpp>

#include "pack/inter_particle_packer.h"
#include "inter_atom_list.h"
#include "../utils/mpi_domain.h"
#include "../utils/mpi_data_types.h"
#include "lattice/ws_utils.h"

InterAtomList::InterAtomList() : nlocalinter(0), nghostinter(0){}

void InterAtomList::makeIndex(AtomList *atom_list, const comm::Domain *p_domain) {
    inter_map.clear();
    // index inter atoms
    for (_type_inter_list::iterator inter_it = inter_list.begin(); inter_it != inter_list.end(); inter_it++) {
        _type_atom_index coords[DIMENSION]; // coords in box (starting from ghost area.)
        ws::getNearLatCoord(*inter_it, p_domain, coords);
        // todo checkout coords boundary.
        const _type_atom_index _atom_near_index = atom_list->lattice.IndexOf3DIndex(coords[0], coords[1], coords[2]);
        inter_map.insert(std::make_pair(_atom_near_index, &(*inter_it))); // save address
    }
    // index ghost inter atoms
    for (_type_inter_list::iterator inter_it = inter_ghost_list.begin();
         inter_it != inter_ghost_list.end(); inter_it++) {
        _type_atom_index coords[DIMENSION]; // coords in box (starting from ghost area.)
        ws::getNearLatCoord(*inter_it, p_domain, coords);
        const _type_atom_index _atom_near_index = atom_list->lattice.IndexOf3DIndex(coords[0], coords[1], coords[2]);
        inter_map.insert(std::make_pair(_atom_near_index, &(*inter_it))); // save address
    }
}

void InterAtomList::addInterAtom(AtomElement atom) {
    inter_list.emplace_back(atom);
    nlocalinter++;
}

AtomElement *InterAtomList::addGhostAtom(AtomElement &ghost_atom) {
    inter_ghost_list.push_back(ghost_atom);
    nghostinter++;
    return &(inter_ghost_list.back());
}

_type_inter_list::iterator InterAtomList::removeInter(_type_inter_list::iterator inter_it) {
    // remove it.
    nlocalinter--;
    return inter_list.erase(inter_it);
}

void InterAtomList::clearGhost() {
    inter_ghost_list.clear(); // clear ghost inter atoms.
    nghostinter = 0; // todo nghostinter is not used.
}
