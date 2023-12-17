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

/**
 * calculate coordinate of nearest lattice for a atom:
 * X,Y,Z is the position of atom and LC is lattice const.
 */
#define VORONOY(X, Y, Z, LC) \
/* first, we get coordinate of major lattice. */ \
_type_atom_index lat_coord_x = static_cast<_type_atom_index>(lround(X / LC)); \
_type_atom_index lat_coord_y = static_cast<_type_atom_index>(lround(Y / LC)); \
_type_atom_index lat_coord_z = static_cast<_type_atom_index>(lround(Z / LC)); \
\
/* second, get delta value which will belongs [-0.5, 0.5]. */ \
const double delta_x = X / LC - lat_coord_x; \
const double delta_y = Y / LC - lat_coord_y; \
const double delta_z = Z / LC - lat_coord_z; \
\
const unsigned int flag_ = (delta_z > 0 ? 4u : 0u) | (delta_y > 0 ? 2u : 0u) | (delta_x > 0 ? 1u : 0u); \
\
/* third, shortest distance calculation. */  \
/* note: we reuse the variable lat_coord_x, lat_coord_y and lat_coord_z as return value.  */ \
lat_coord_x = 2 * lat_coord_x;                  \
if (normal_vector[flag_][0] * delta_x +          \
    normal_vector[flag_][1] * delta_y +              \
    normal_vector[flag_][2] * delta_z + d >= 0.0) {  \
    /* belongs to major lattice. */ \
    lat_coord_x += offset[flag_][0]; \
    lat_coord_y += offset[flag_][1]; \
    lat_coord_z += offset[flag_][2]; \
}

const box::_type_flag_32 ws::isOutBox(_type_atom_location src_atom_x[DIMENSION], const comm::Domain *p_domain) {
    VORONOY(src_atom_x[0], src_atom_x[1], src_atom_x[2], p_domain->lattice_const)

    lat_coord_x -= 2 * p_domain->sub_box_lattice_region.x_low;
    lat_coord_y -= p_domain->sub_box_lattice_region.y_low;
    lat_coord_z -= p_domain->sub_box_lattice_region.z_low;

    box::_type_flag_32 flag = box::IN_BOX;
    if (lat_coord_x < 0) {
        flag |= box::OUT_BOX_X_LITTER;
    } else if (lat_coord_x >= 2 * p_domain->sub_box_lattice_size[0]) {
        flag |= box::OUT_BOX_X_BIG;
    }
    if (lat_coord_y < 0) {
        flag |= box::OUT_BOX_Y_LITTER;
    } else if (lat_coord_y >= p_domain->sub_box_lattice_size[1]) {
        flag |= box::OUT_BOX_Y_BIG;
    }
    if (lat_coord_z < 0) {
        flag |= box::OUT_BOX_Z_LITTER;
    } else if (lat_coord_z >= p_domain->sub_box_lattice_size[2]) {
        flag |= box::OUT_BOX_Z_BIG;
    }
    return flag;
}

_type_atom_index ws::findNearLatAtom(AtomList *atom_list, const AtomElement &src_atom, const comm::Domain *p_domain) {
    _type_atom_index coords[DIMENSION];
    getNearLatCoord(src_atom, p_domain, coords);
    _type_atom_index near_index = atom_list->lattice.IndexOf3DIndex(coords[0], coords[1], coords[2]);
    return near_index;
}

_type_atom_index ws::findNearLatAtomInSubBox(AtomList *atom_list, const _type_atom_location src_x[DIMENSION],
                                             const comm::Domain *p_domain) {
    _type_atom_index near_index = findNearLatIndexInSubBox(atom_list->lattice, src_x, p_domain);
    return near_index;
}

_type_atom_index ws::findNearLatIndexInSubBox(const BccLattice &lattice, const _type_atom_location src_x[DIMENSION],
                                              const comm::Domain *p_domain) {
    // get the lattice coordinate of nearest atom first.
    _type_atom_index coord_to_sub_box[DIMENSION];
    getNearLatSubBoxCoord(src_x, p_domain, coord_to_sub_box);
    _type_atom_index &j = coord_to_sub_box[0];
    _type_atom_index &k = coord_to_sub_box[1];
    _type_atom_index &l = coord_to_sub_box[2];

    // if nearest atom is out of box.
    if (j < 0 || k < 0 || l < 0 ||
        j >= 2 * p_domain->sub_box_lattice_size[0] ||
        k >= p_domain->sub_box_lattice_size[1] ||
        l >= p_domain->sub_box_lattice_size[2]) {
        return box::IndexNotExists;
    }
    // calculate atom index in ghost included sub-box
    j += 2 * p_domain->lattice_size_ghost[0];
    k += p_domain->lattice_size_ghost[1];
    l += p_domain->lattice_size_ghost[2];
    return lattice.IndexOf3DIndex(j, k, l);
}

void ws::getNearLatCoord(const AtomElement &src_atom, const comm::Domain *p_domain,
                         _type_atom_index coords[DIMENSION]) {
    VORONOY(src_atom.x[0], src_atom.x[1], src_atom.x[2], p_domain->lattice_const)

    coords[0] = lat_coord_x - 2 * p_domain->ghost_ext_lattice_region.x_low;
    coords[1] = lat_coord_y - p_domain->ghost_ext_lattice_region.y_low;
    coords[2] = lat_coord_z - p_domain->ghost_ext_lattice_region.z_low;
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
 * otherwise return lattice coordinate of minor lattice.
 */
void ws::getNearLatSubBoxCoord(const _type_atom_location src_x[DIMENSION], const comm::Domain *p_domain,
                               _type_atom_index coords[DIMENSION]) {
    VORONOY(src_x[0], src_x[1], src_x[2], p_domain->lattice_const)

    // get relative coords
    coords[0] = lat_coord_x - 2 * p_domain->sub_box_lattice_region.x_low;
    coords[1] = lat_coord_y - p_domain->sub_box_lattice_region.y_low;
    coords[2] = lat_coord_z - p_domain->sub_box_lattice_region.z_low;
}

bool ws::isInBox(const AtomElement &src_atom, const comm::Domain *p_domain) {
    return isInBox(src_atom.x[0], src_atom.x[1], src_atom.x[2], p_domain);
}

bool ws::isInBox(const double rx, const double ry, const double rz, const comm::Domain *p_domain) {
    VORONOY(rx, ry, rz, p_domain->lattice_const)

    lat_coord_x -= 2 * p_domain->sub_box_lattice_region.x_low;
    lat_coord_y -= p_domain->sub_box_lattice_region.y_low;
    lat_coord_z -= p_domain->sub_box_lattice_region.z_low;
    return (lat_coord_x >= 0 && lat_coord_y >= 0 && lat_coord_z >= 0 &&
            lat_coord_x < 2 * p_domain->sub_box_lattice_size[0] &&
            lat_coord_y < p_domain->sub_box_lattice_size[1] &&
            lat_coord_z < p_domain->sub_box_lattice_size[2]);
}
