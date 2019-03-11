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

AtomElement &ws::findNearLatAtom(AtomList *atom_list, const AtomElement &src_atom, const Domain *p_domain) {
    _type_atom_index coords[DIMENSION];
    getNearLatCoord(src_atom, p_domain, coords);
    _type_atom_index near_index = atom_list->IndexOf3DIndex(coords[0], coords[1], coords[2]);
    return atom_list->getAtomEleByLinearIndex(near_index); // todo return _atoms[l][k][j];
}

AtomElement *ws::findNearLatAtomInSubBox(AtomList *atom_list, const AtomElement &src_atom, const Domain *p_domain) {
    _type_atom_index near_index = findNearLatIndexInSubBox(atom_list, src_atom, p_domain);
    if (near_index == box::IndexNotExists) {
        return nullptr;
    }
    return &(atom_list->getAtomEleByLinearIndex(near_index));
}

_type_atom_index
ws::findNearLatIndexInSubBox(AtomList *atom_list, const AtomElement &src_atom, const Domain *p_domain) {
    // get the lattice coordinate of nearest atom first.
    _type_atom_index coord_to_sub_box[DIMENSION];
    getNearLatSubBoxCoord(src_atom, p_domain, coord_to_sub_box);
    _type_atom_index &j = coord_to_sub_box[0];
    _type_atom_index &k = coord_to_sub_box[1];
    _type_atom_index &l = coord_to_sub_box[2];

    // if nearest atom is out of box.
    if (j < 0 || k < 0 || l < 0 ||
        j >= 2 * p_domain->lattice_size_sub_box[0] ||
        k >= p_domain->lattice_size_sub_box[1] ||
        l >= p_domain->lattice_size_sub_box[2]) {
        return box::IndexNotExists;
    }
    // calculate atom index in ghost included sub-box
    j += 2 * p_domain->lattice_size_ghost[0];
    k += p_domain->lattice_size_ghost[1];
    l += p_domain->lattice_size_ghost[2];
    return atom_list->IndexOf3DIndex(j, k, l);
}

void ws::getNearLatCoord(const AtomElement &src_atom, const Domain *p_domain, _type_atom_index coords[DIMENSION]) {
    auto j = static_cast<_type_atom_index>(lround(src_atom.x[0] * 2 / p_domain->lattice_const));
    auto k = static_cast<_type_atom_index>(lround(src_atom.x[1] * 2 / p_domain->lattice_const));
    auto l = static_cast<_type_atom_index>(lround(src_atom.x[2] * 2 / p_domain->lattice_const));
    k = k < 0 ? k / 2 - 1 : k / 2;
    l = l < 0 ? l / 2 - 1 : l / 2;
    coords[0] = j - 2 * p_domain->lattice_coord_ghost_region.x_low;
    coords[1] = k - p_domain->lattice_coord_ghost_region.y_low;
    coords[2] = l - p_domain->lattice_coord_ghost_region.z_low;
}

void ws::getNearLatSubBoxCoord(const AtomElement &src_atom, const Domain *p_domain,
                               _type_atom_index coords[DIMENSION]) {
    auto j = static_cast<_type_atom_index>(lround(src_atom.x[0] * 2 / p_domain->lattice_const));
    auto k = static_cast<_type_atom_index>(lround(src_atom.x[1] * 2 / p_domain->lattice_const));
    auto l = static_cast<_type_atom_index>(lround(src_atom.x[2] * 2 / p_domain->lattice_const));
    k = k < 0 ? k / 2 - 1 : k / 2;
    l = l < 0 ? l / 2 - 1 : l / 2;
    coords[0] = j - 2 * p_domain->lattice_coord_sub_box_region.x_low;
    coords[1] = k - p_domain->lattice_coord_sub_box_region.y_low;
    coords[2] = l - p_domain->lattice_coord_sub_box_region.z_low;
}

bool ws::isInBox(const AtomElement &src_atom, const Domain *p_domain) {
    auto j = static_cast<_type_atom_index>(lround(src_atom.x[0] * 2 / p_domain->lattice_const));
    auto k = static_cast<_type_atom_index>(lround(src_atom.x[1] * 2 / p_domain->lattice_const));
    auto l = static_cast<_type_atom_index>(lround(src_atom.x[2] * 2 / p_domain->lattice_const));
    k = k / 2;
    l = l / 2;
    j -= 2 * p_domain->lattice_coord_sub_box_region.x_low;
    k -= p_domain->lattice_coord_sub_box_region.y_low;
    l -= p_domain->lattice_coord_sub_box_region.z_low;
    return (j >= 0 && k >= 0 && l >= 0 &&
            j < p_domain->lattice_size_sub_box[0] &&
            k < p_domain->lattice_size_sub_box[1] &&
            l < p_domain->lattice_size_sub_box[2]);
}
