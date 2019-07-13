//
// Created by genshen on 2019-01-02.
//

#include <cmath>
#include "ws_utils.h"

// the consts are only used in this cpp file.
namespace ws {
    /**
     * The coordinate offset of 8 minor lattices surrounding a major lattice.
     * If the major lattice is determined, we can get the lattice coordinate of a minor lattice
     * by this offset.
     */
    const _type_atom_index offset[8][DIMENSION] = { //offset index: 0-> x, 1 -> y, 2 -> z
            {-1, -1, -1}, // x low,  y low,  z low, (b000)
            {1,  -1, -1}, // x high, y low,  z low, (b001)
            {-1, 0,  -1}, // x low,  y high, z low, (b010)
            {1,  0,  -1}, // x high, y high, z low, (b011)
            {-1, -1, 0},  // x low,  y low,  z high,(b100)
            {1,  -1, 0},  // x high, y low,  z high,(b101)
            {-1, 0,  0},  // x low,  y high, z high,(b110)
            {1,  0,  0},  // x high, y high, z high,(b111)
    };

    /**
     * Assuming plane: Ax+By+Cz+d = 0.
     * For major lattice with coordinate: (0,0,0), it has 8 1nn lattices.
     * Then, we can compute the 8 vertical bisector planes between major lattice and each 1nn lattices:
     * 0: -x-y-z-3/4 = 0;
     * 1: x-y-z-3/4 = 0;
     * 2: -x+y-z-3/4 = 0;
     * 3: x+y-z-3/4 = 0;
     * 4: -x-y+z-3/4 = 0;
     * 5: x-y+z-3/4 = 0;
     * 6: -x+y+z-3/4 = 0;
     * 7: x+y+z-3/4 = 0;
     *
     * This array saves the normal vectors of 8 vertical bisector planes.
     *
     * Noticing that, for major lattice position:(0, 0, 0), we always have A*0+B*0+C*0+d < 0.
     */
    const double normal_vector[8][3] = {
            {-1, -1, -1}, // normal vector of vertical bisector plane of (-0.5, -0.5, -0.5) and (0,0,0)
            {1,  -1, -1}, // (0.5,  -0.5, -0.5) and (0,0,0)
            {-1, 1,  -1}, // (-0.5, 0.5, -0.5)  and (0,0,0)
            {1,  1,  -1}, // (0.5, 0.5, -0.5)  and (0,0,0)
            {-1, -1, 1},  // (-0.5, -0.5, 0.5)  and (0,0,0)
            {1,  -1, 1},  // (0.5,  -0.5, 0.5)  and (0,0,0)
            {-1, 1,  1},  // (-0.5, 0.5,  0.5)  and (0,0,0)
            {1,  1,  1},  // (0.5,  0.5, 0.5)  and (0,0,0)
    };

    /**
     * d = -(A*x_m + B*y_m + C*z_m), in which (A,B,C) is normal vector (such as (-1, -1, -1) ),
     * (x_m, y_m, z_m) is middle pointer (such as (-1/4, -1/4, -1/4)).
     */
    const double d = -3.0 / 4.0;

}

const box::_type_flag_32 ws::isOutBox(const AtomElement &src_atom, const comm::Domain *p_domain) {
    auto j = static_cast<_type_atom_index>(lround(src_atom.x[0] * 2 / p_domain->lattice_const));
    auto k = static_cast<_type_atom_index>(lround(
            (src_atom.x[1] - (j % 2) * p_domain->lattice_const / 2) / p_domain->lattice_const));
    auto l = static_cast<_type_atom_index>(lround(
            (src_atom.x[2] - (j % 2) * p_domain->lattice_const / 2) / p_domain->lattice_const));
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

AtomElement &ws::findNearLatAtom(AtomList *atom_list, const AtomElement &src_atom, const comm::Domain *p_domain) {
    _type_atom_index coords[DIMENSION];
    getNearLatCoord(src_atom, p_domain, coords);
    _type_atom_index near_index = atom_list->lattice.IndexOf3DIndex(coords[0], coords[1], coords[2]);
    return atom_list->getAtomEleByLinearIndex(near_index); // todo return _atoms[l][k][j];
}

AtomElement *ws::findNearLatAtomInSubBox(AtomList *atom_list, const AtomElement &src_atom,
                                         const comm::Domain *p_domain) {
    _type_atom_index near_index = findNearLatIndexInSubBox(atom_list, src_atom, p_domain);
    if (near_index == box::IndexNotExists) {
        return nullptr;
    }
    return &(atom_list->getAtomEleByLinearIndex(near_index));
}

_type_atom_index ws::findNearLatIndexInSubBox(AtomList *atom_list, const AtomElement &src_atom,
                                              const comm::Domain *p_domain) {
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
    return atom_list->lattice.IndexOf3DIndex(j, k, l);
}

void ws::getNearLatCoord(const AtomElement &src_atom, const comm::Domain *p_domain,
                         _type_atom_index coords[DIMENSION]) {
    auto j = static_cast<_type_atom_index>(lround(src_atom.x[0] * 2 / p_domain->lattice_const));
    auto k = static_cast<_type_atom_index>(lround(
            (src_atom.x[1] - (j % 2) * p_domain->lattice_const / 2) / p_domain->lattice_const));
    auto l = static_cast<_type_atom_index>(lround(
            (src_atom.x[2] - (j % 2) * p_domain->lattice_const / 2) / p_domain->lattice_const));

    coords[0] = j - 2 * p_domain->lattice_coord_ghost_region.x_low;
    coords[1] = k - p_domain->lattice_coord_ghost_region.y_low;
    coords[2] = l - p_domain->lattice_coord_ghost_region.z_low;
}

/**
 * for example, for a position x = (5.6, 7.01, 8.9) and lattice const is c = 1.0;
 * 1. we can find coordinate of major lattice M = (6, 7, 9) by: (int) round(x/c).
 * so, the belonging lattice can be in: M  and 8 1nn neighbors lattices of M.
 *
 * 2. determine minor lattice by delta distance of atom position and major lattice position.
 * if delta_x < 0 ? minor lattice at back face, otherwise minor lattice at front face;
 * if delta_y < 0 ? minor lattice at left face, otherwise minor lattice at right face;
 * if delta_z < 0 ? minor lattice at bottom face, otherwise minor lattice at top face;
 *
 * 3. shortest distance calculation:
 * calculate distance (denoted as d1) of atom position to major lattice and
 * distance (denoted as d2) of atom position to minor lattice.
 * If d1 < d2, return lattice coordinate of major lattice,
 * otherwise return lattice coordinate of miinor lattice.
 */
void ws::getNearLatSubBoxCoord(const AtomElement &src_atom, const comm::Domain *p_domain,
                               _type_atom_index coords[DIMENSION]) {
    // first, we get coordinate of major lattice.
    _type_atom_index major_lat_x = static_cast<_type_atom_index>(lround(src_atom.x[0] / p_domain->lattice_const));
    _type_atom_index major_lat_y = static_cast<_type_atom_index>(lround(src_atom.x[1] / p_domain->lattice_const));
    _type_atom_index major_lat_z = static_cast<_type_atom_index>(lround(src_atom.x[2] / p_domain->lattice_const));

    // second, get delta value which will belongs [-0.5, 0.5].
    const double delta_x = src_atom.x[0] / p_domain->lattice_const - major_lat_x;
    const double delta_y = src_atom.x[1] / p_domain->lattice_const - major_lat_y;
    const double delta_z = src_atom.x[2] / p_domain->lattice_const - major_lat_z;

    const unsigned int flag = (delta_z > 0 ? 4u : 0u) | (delta_y > 0 ? 2u : 0u) | (delta_x > 0 ? 1u : 0u);

    // third, shortest distance calculation.
    if (normal_vector[flag][0] * delta_x +
        normal_vector[flag][1] * delta_y +
        normal_vector[flag][2] * delta_z + d < 0.0) {
        // belongs to major lattice
        coords[0] = 2 * major_lat_x;
        coords[1] = major_lat_y;
        coords[2] = major_lat_z;
    } else {
        coords[0] = 2 * major_lat_x + offset[flag][0];
        coords[1] = major_lat_y + offset[flag][1];
        coords[2] = major_lat_z + offset[flag][2];
    }
    // relative coords
    coords[0] -= 2 * p_domain->lattice_coord_sub_box_region.x_low;
    coords[1] -= p_domain->lattice_coord_sub_box_region.y_low;
    coords[2] -= p_domain->lattice_coord_sub_box_region.z_low;
}

bool ws::isInBox(const AtomElement &src_atom, const comm::Domain *p_domain) {
    auto j = static_cast<_type_atom_index>(lround(src_atom.x[0] * 2 / p_domain->lattice_const));
    auto k = static_cast<_type_atom_index>(lround(
            (src_atom.x[1] - (j % 2) * p_domain->lattice_const / 2) / p_domain->lattice_const));
    auto l = static_cast<_type_atom_index>(lround(
            (src_atom.x[2] - (j % 2) * p_domain->lattice_const / 2) / p_domain->lattice_const));

    j -= 2 * p_domain->lattice_coord_sub_box_region.x_low;
    k -= p_domain->lattice_coord_sub_box_region.y_low;
    l -= p_domain->lattice_coord_sub_box_region.z_low;
    return (j >= 0 && k >= 0 && l >= 0 &&
            j < p_domain->lattice_size_sub_box[0] &&
            k < p_domain->lattice_size_sub_box[1] &&
            l < p_domain->lattice_size_sub_box[2]);
}
