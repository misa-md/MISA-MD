//
// Created by genshen on 2018-12-31.
//

#include <types/pre_define.h>
#include <cmath>
#include <hardware_accelerate.hpp>
#include "atom_set.h"

AtomSet::AtomSet(Domain *domain, double latticeconst, double cutoffRadiusFactor)
        : p_domain(domain), _latticeconst(latticeconst),
          _cutoffRadius(cutoffRadiusFactor * latticeconst){

    _cutlattice = static_cast<int>(ceil(cutoffRadiusFactor));

    numberoflattice = p_domain->lattice_size_ghost_extended[0] * p_domain->lattice_size_ghost_extended[1] *
                      p_domain->lattice_size_ghost_extended[2];
    // printf("number:%d, %d, %d, %d, %d\n", numberoflattice,
    // p_domain->getSubBoxLatticeSize(0), p_domain->getSubBoxLatticeSize(1), p_domain->getSubBoxLatticeSize(2),
    // p_domain->getSubBoxLatticeSize(0)*p_domain->getSubBoxLatticeSize(1)*p_domain->getSubBoxLatticeSize(2));
    atom_list = new AtomList(p_domain->lattice_size_ghost_extended[0],
                             p_domain->lattice_size_ghost_extended[1],
                             p_domain->lattice_size_ghost_extended[2],
                             p_domain->lattice_size_sub_box[0],
                             p_domain->lattice_size_sub_box[1],
                             p_domain->lattice_size_sub_box[2],
                             p_domain->lattice_size_ghost[0],
                             p_domain->lattice_size_ghost[1],
                             p_domain->lattice_size_ghost[2]);

    inter_atom_list = new InterAtomList();

    if (isAccelerateSupport()) {
        accelerateInit(p_domain->lattice_coord_sub_box_region.x_low,
                       p_domain->lattice_coord_sub_box_region.y_low,
                       p_domain->lattice_coord_sub_box_region.z_low,
                       p_domain->lattice_size_sub_box[0],
                       p_domain->lattice_size_sub_box[1],
                       p_domain->lattice_size_sub_box[2],
                       p_domain->lattice_coord_ghost_region.x_low,
                       p_domain->lattice_coord_ghost_region.y_low,
                       p_domain->lattice_coord_ghost_region.z_low,
                       p_domain->lattice_size_ghost_extended[0],
                       p_domain->lattice_size_ghost_extended[1],
                       p_domain->lattice_size_ghost_extended[2]);
    }
}

AtomSet::~AtomSet() {
    delete atom_list;
    delete inter_atom_list;
}

void AtomSet::calculateNeighbourIndices() {
    double x, y, z;
    int mark = 0;
    double cut_times_lattice = _cutoffRadius / _latticeconst; // todo use cutoffRadiusFactor.
    std::vector<_type_atom_index>::iterator neighbourOffsetsIter;
    for (int zIndex = -_cutlattice; zIndex <= _cutlattice; zIndex++) { // loop for (2*_cutlattice + 1) times.
        for (int yIndex = -_cutlattice; yIndex <= _cutlattice; yIndex++) {
            for (int xIndex = -_cutlattice * 2; xIndex <= _cutlattice * 2; xIndex++) {
                // 体心
                z = (double) zIndex + (((double) (xIndex % 2)) / 2); // zIndex plus 1/2 (odd) or 0(even).
                y = (double) yIndex + (((double) (xIndex % 2)) / 2);
                x = (double) xIndex / 2;
                _type_atom_index offset;
                double r = sqrt(x * x + y * y + z * z);
                if (r < (cut_times_lattice + 0.4)) { // todo 0.4?
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
                if (r < (cut_times_lattice + 0.4)) {
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

void AtomSet::addAtom(_type_atom_id id, double rx, double ry, double rz, double vx, double vy, double vz) {
    int i;
    if ((rx >= p_domain->meas_sub_box_region.x_low) &&
        (rx < p_domain->meas_sub_box_region.x_high) &&
        (ry >= p_domain->meas_sub_box_region.y_low) &&
        (ry < p_domain->meas_sub_box_region.y_high) &&
        (rz >= p_domain->meas_sub_box_region.z_low) &&
        (rz < p_domain->meas_sub_box_region.z_high)) {
        int lattice[3];
        lattice[0] = rx * 2 / _latticeconst + 0.5;
        lattice[1] = ry * 2 / _latticeconst + 0.5;
        lattice[2] = rz * 2 / _latticeconst + 0.5;
        lattice[1] = lattice[1] / 2;
        lattice[2] = lattice[2] / 2;
        lattice[0] -= p_domain->lattice_coord_ghost_region.x_low;
        lattice[1] -= p_domain->lattice_coord_ghost_region.y_low;
        lattice[2] -= p_domain->lattice_coord_ghost_region.z_low;
        i = (((p_domain->lattice_size_ghost_extended[1])) * lattice[2] + lattice[1]) *
            ((p_domain->lattice_size_ghost_extended[0])) + lattice[0];
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

_type_atom_count AtomSet::getnlocalatom() {
    return (p_domain->lattice_size_sub_box[0] * p_domain->lattice_size_sub_box[1] *
            p_domain->lattice_size_sub_box[2]);
}
