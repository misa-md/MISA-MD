//
// Created by genshen on 2019-01-02.
//

#include <cmath>
#include "ws_utils.h"

const box::_type_flag_32 ws::isOutBox(const AtomElement &src_atom, const Domain *p_domain) {
    auto j = static_cast<_type_atom_index>(lround(src_atom.x[0] * 2 / p_domain->lattice_const));
    auto k = static_cast<_type_atom_index>(lround(src_atom.x[1] * 2 / p_domain->lattice_const));
    auto l = static_cast<_type_atom_index>(lround(src_atom.x[2] * 2 / p_domain->lattice_const));
    k = k < 0 ? k / 2 - 1 : k / 2; // for example: right k=-1, then -1/2 = 0, but we expect left k=-1
    l = l < 0 ? l / 2 - 1 : l / 2;
    j -= 2 * p_domain->lattice_coord_sub_box_region.x_low;
    k -= p_domain->lattice_coord_sub_box_region.y_low;
    l -= p_domain->lattice_coord_sub_box_region.z_low;

    box::_type_flag_32 flag = box::IN_BOX;
    if (j < 0) {
        flag |= box::OUT_BOX_X_LITTER;
    } else if (j >= 2 * p_domain->lattice_size_sub_box[0]) {
        flag |= box::OUT_BOX_X_BIG;
    }
    if (k < 0) {
        flag |= box::OUT_BOX_Y_LITTER;
    } else if (k >= p_domain->lattice_size_sub_box[1]) {
        flag |= box::OUT_BOX_Y_BIG;
    }
    if (l < 0) {
        flag |= box::OUT_BOX_Z_LITTER;
    } else if (l >= p_domain->lattice_size_sub_box[2]) {
        flag |= box::OUT_BOX_Z_BIG;
    }
    return flag;
}
