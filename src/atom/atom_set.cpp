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
        : _cutoffRadius(cutoff_radius){
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
}

AtomSet::~AtomSet() {
    delete atom_list;
    delete inter_atom_list;
}

void AtomSet::calcNeighbourIndices(const double cutoff_radius_factor, const _type_lattice_size cut_lattice) {
    double x, y, z;
    int mark = 0;
    std::vector<_type_atom_index>::iterator neighbourOffsetsIter;
    for (_type_atom_index zIndex = -cut_lattice;
         zIndex <= cut_lattice; zIndex++) { // loop for (2*_cutlattice + 1) times.
        for (_type_atom_index yIndex = -cut_lattice; yIndex <= cut_lattice; yIndex++) {
            for (_type_atom_index xIndex = -cut_lattice * 2; xIndex <= cut_lattice * 2; xIndex++) {
                // 体心
                z = (double) zIndex + (((double) (xIndex % 2)) / 2); // zIndex plus 1/2 (odd) or 0(even).
                y = (double) yIndex + (((double) (xIndex % 2)) / 2);
                x = (double) xIndex / 2;
                _type_atom_index offset;
                double r = sqrt(x * x + y * y + z * z);
                if (r < (cutoff_radius_factor + 0.4)) { // todo 0.4?
                    offset = atom_list->IndexOf3DIndex(xIndex, yIndex, zIndex);
                    if (offset > 0) {
                        NeighbourOffsets.push_back(offset);
                    }
                }

                // 晶格点
                z = (double) zIndex - (((double) (xIndex % 2)) / 2);
                y = (double) yIndex - (((double) (xIndex % 2)) / 2);
                x = (double) xIndex / 2;
                r = sqrt(x * x + y * y + z * z);
                if (r < (cutoff_radius_factor + 0.4)) {
                    offset = atom_list->IndexOf3DIndex(xIndex, yIndex, zIndex);
                    if (offset > 0) {
                        for (neighbourOffsetsIter = NeighbourOffsets.begin();
                             neighbourOffsetsIter != NeighbourOffsets.end(); neighbourOffsetsIter++) {
                            if (*neighbourOffsetsIter == offset) {
                                mark = 1;
                            }
                        }
                        if (mark != 1) {
                            NeighbourOffsets.push_back(offset);
                        }
                        mark = 0;
                    }
                }
            }
        }
    }
}

void
AtomSet::addAtom(comm::Domain *p_domain, _type_atom_id id,
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
