//
// Created by genshen on 2018-12-31.
//

#include <cmath>
#include <types/pre_define.h>
#include <hardware_accelerate.hpp>
#include <domain/domain.h>

#include "atom_set.h"

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
    neighbours = new NeighbourIndex<AtomElement>(*atom_list);
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
AtomSet::addAtom(comm::BccDomain *p_domain, _type_atom_id id,
                 double rx, double ry, double rz, double vx, double vy, double vz) {
    int i;
    if ((rx >= p_domain->meas_sub_box_region.x_low) &&
        (rx < p_domain->meas_sub_box_region.x_high) &&
        (ry >= p_domain->meas_sub_box_region.y_low) &&
        (ry < p_domain->meas_sub_box_region.y_high) &&
        (rz >= p_domain->meas_sub_box_region.z_low) &&
        (rz < p_domain->meas_sub_box_region.z_high)) {
        int lattice[3];
        lattice[0] = rx * 2 / p_domain->lattice_const + 0.5;
        lattice[1] = ry * 2 / p_domain->lattice_const + 0.5;
        lattice[2] = rz * 2 / p_domain->lattice_const + 0.5;
        lattice[1] = lattice[1] / 2;
        lattice[2] = lattice[2] / 2;
        lattice[0] -= p_domain->dbx_lattice_coord_ghost_region.x_low;
        lattice[1] -= p_domain->dbx_lattice_coord_ghost_region.y_low;
        lattice[2] -= p_domain->dbx_lattice_coord_ghost_region.z_low;
        i = (((p_domain->dbx_lattice_size_ghost_extended[1])) * lattice[2] + lattice[1]) *
            ((p_domain->dbx_lattice_size_ghost_extended[0])) + lattice[0];
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(i);
        atom_.id = id;
        atom_.x[0] = rx;
        atom_.x[1] = ry;
        atom_.x[2] = rz;
        atom_.v[0] = vx;
        atom_.v[1] = vy;
        atom_.v[2] = vz;
    }
}

_type_atom_count AtomSet::getnlocalatom(comm::Domain *p_domain) {
    return (p_domain->lattice_size_sub_box[0] * p_domain->lattice_size_sub_box[1] *
            p_domain->lattice_size_sub_box[2]);
}

#ifdef MD_DEV_MODE

_type_atom_count AtomSet::realAtoms() {
    _type_atom_count count = 0;
    atom_list->foreachSubBoxAtom(
            [&count](AtomElement &_atom_ref) {
                if (_atom_ref.type != atom_type::INVALID) {
                    count++;
                }
            }
    );
    return count;
}

#endif
