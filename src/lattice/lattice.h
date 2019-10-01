//
// Created by genshen on 2019-04-20.
//

#ifndef CRYSTALMD_LATTICE_H
#define CRYSTALMD_LATTICE_H


#include "types/pre_define.h"

class BccLattice {
public:
    // the total lattice count.
    const _type_atom_count _size;
    // note: _size_x is twice times than the lattice size in x dimension.
    const _type_atom_count _size_x, _size_y, _size_z;
    const _type_atom_count _size_sub_box_x, _size_sub_box_y, _size_sub_box_z;
    const _type_atom_count purge_ghost_count_x, purge_ghost_count_y, purge_ghost_count_z;

public:
    /**
     *
     * initialize atom lattice with the size of atoms(including ghost atoms).
     * @param size_x,size_y,size_z the atoms count each dimension(including ghost atoms).
     * @param ghost_count_x,ghost_count_y, ghost_count_z the count of ghost atoms at left or right of each dimension.
     * (in fact left and right has same count of ghost atoms)
     * @note the size at x dimension is doubled (including @param size_x,ghost_count_x).
     */
    BccLattice(const _type_atom_count size_x,
               const _type_atom_count size_y,
               const _type_atom_count size_z,
               const _type_atom_count size_sub_box_x,
               const _type_atom_count size_sub_box_y,
               const _type_atom_count size_sub_box_z,
               const _type_atom_count ghost_count_x,
               const _type_atom_count ghost_count_y,
               const _type_atom_count ghost_count_z);

    /**
     * get the atom  3d index by linear index.
     * index = (zIndex * p_domain->getGhostLatticeSize(1) + yIndex) * p_domain->getGhostLatticeSize(0) + xIndex;
     * @param index target linear index
     * @param x 3d index in x dimension
     * @param y 3d index in y dimension
     * @param z 3d index in z dimension
     */
    inline void get3DIndexByLinearIndex(_type_atom_index index, _type_atom_index &x,
                                        _type_atom_index &y, _type_atom_index &z) const {
        x = index % _size_x;
        index = index / _size_x;
        y = index % _size_y;
        z = index / _size_y;
    }

    /**
     * get linear index of 3d atoms array
     * @param xIndex index at x dimension
     * @param yIndex index at y dimension
     * @param zIndex
     * @return
     */
    inline _type_atom_index
    IndexOf3DIndex(_type_atom_index xIndex, _type_atom_index yIndex, _type_atom_index zIndex) const {
        return (zIndex * _size_y + yIndex) * _size_x + xIndex;
    }
};


#endif //CRYSTALMD_LATTICE_H
