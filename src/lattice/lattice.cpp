//
// Created by genshen on 2019-04-20.
//

#include "lattice.h"

BccLattice::BccLattice(const _type_atom_count size_x, const _type_atom_count size_y, const _type_atom_count size_z,
                       const _type_atom_count size_sub_box_x, const _type_atom_count size_sub_box_y,
                       const _type_atom_count size_sub_box_z, const _type_atom_count ghost_count_x,
                       const _type_atom_count ghost_count_y, const _type_atom_count ghost_count_z)
        : _size(size_x * size_y * size_z),
          _size_x(size_x), _size_y(size_y), _size_z(size_z),
          _size_sub_box_x(size_sub_box_x),
          _size_sub_box_y(size_sub_box_y),
          _size_sub_box_z(size_sub_box_z),
          purge_ghost_count_x(ghost_count_x),
          purge_ghost_count_y(ghost_count_y),
          purge_ghost_count_z(ghost_count_z) {

}
