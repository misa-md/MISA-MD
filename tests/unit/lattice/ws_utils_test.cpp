//
// Created by genshen on 2019-04-21.
//

#include <gtest/gtest.h>
#include <lattice/ws_utils.h>


/**
 * tests for BCC lattice ws(Wigner-Seitz) method.
 */

TEST(ws_utils_getNearLatSubBoxCoord_test, _ws_utils_test) {
    const int grid_size[3] = {2, 2, 2};
    const int grid_coord[3] = {0, 0, 0};
    const int64_t space[3] = {50 * grid_size[0], 60 * grid_size[1], 72 * grid_size[2]};
    const double lattice_const = 0.86;
    const double cutoff_radius_factor = 1.1421;
    comm::Domain *p_domain = comm::Domain::Builder()
            .setPhaseSpace(space)
            .setCutoffRadius(cutoff_radius_factor)
            .setLatticeConst(lattice_const)
            .localBuild(grid_size, grid_coord);

    for (unsigned int z = 0; z < 72; z++) {
        for (unsigned int y = 0; y < 60; y++) {
            for (unsigned int x = 0; x < 2 * 50; x++) {
                AtomElement src_atom;
                src_atom.x[0] = x * (lattice_const / 2);
                src_atom.x[1] = x % 2 == 0 ? lattice_const * y : lattice_const * y + lattice_const / 2;
                src_atom.x[2] = x % 2 == 0 ? lattice_const * z : lattice_const * z + lattice_const / 2;

                _type_atom_index coords[DIMENSION];
                ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
                EXPECT_EQ(coords[0], x);
                EXPECT_EQ(coords[1], y);
                EXPECT_EQ(coords[2], z);

                //  add delta values
                src_atom.x[0] = src_atom.x[0] + 0.1 * lattice_const;
                src_atom.x[1] = src_atom.x[1] - 0.2 * lattice_const;
                src_atom.x[2] = src_atom.x[2] + 0.15 * lattice_const;

                ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
                EXPECT_EQ(coords[0], x);
                EXPECT_EQ(coords[1], y);
                EXPECT_EQ(coords[2], z);
            }
        }
    }
}

// test getNearLatSubBoxCoord when x[0] coord is even(e.g x[0] = 0).
TEST(getNearLatSubBoxCoord_minus_x_even_test, _ws_utils_test) {
    const int grid_size[3] = {2, 2, 2};
    const int grid_coord[3] = {0, 0, 0};
    const int64_t space[3] = {25 * grid_size[0], 25 * grid_size[1], 25 * grid_size[2]};
    const double lattice_const = 2.85532;
    const double cutoff_radius_factor = 1.96125;
    comm::Domain *p_domain = comm::Domain::Builder()
            .setPhaseSpace(space)
            .setCutoffRadius(cutoff_radius_factor)
            .setLatticeConst(lattice_const)
            .localBuild(grid_size, grid_coord);

/**
 * figure of atom at z dimension if x is odd.
      z --->
      1     0.5    0    -0.5    1   -1.5    -2    -2.5   -3
      |____________|____________|____________|____________|
index:1            0            1            2            3
z lattice coord:
  z in [0.5, 1.5] -> index z = 1;
  z in [-0.5, 0.5] -> index z = 0;
  z in [-1.5, -0.5] -> index z = -1;
  z in [-2.5, -1.5] -> index z = -2;
  z in [-3.5, -2.5] -> index z = -3;
 */
    // case 1: make z coord < 0.
    AtomElement src_atom;
    src_atom.x[0] = 0;
    src_atom.x[1] = 0;
    src_atom.x[2] = 0 - 0.01;
    _type_atom_index coords[DIMENSION];
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[2], 0);

    // case 2:
    src_atom.x[2] = lattice_const * (-2) - 0.01;
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[2], -2);

    // case 3:
    src_atom.x[2] = lattice_const * (-2.44);
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[2], -2);

    // case 4:
    src_atom.x[2] = lattice_const * (-2.66);
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[2], -3);
}

// test getNearLatSubBoxCoord when x[0] is odd.
TEST(getNearLatSubBoxCoord_minus_x_odd_test, _ws_utils_test) {
    const int grid_size[3] = {2, 2, 2};
    const int grid_coord[3] = {0, 0, 0};
    const int64_t space[3] = {25 * grid_size[0], 25 * grid_size[1], 25 * grid_size[2]};
    const double lattice_const = 2.85532;
    const double cutoff_radius_factor = 1.96125;
    comm::Domain *p_domain = comm::Domain::Builder()
            .setPhaseSpace(space)
            .setCutoffRadius(cutoff_radius_factor)
            .setLatticeConst(lattice_const)
            .localBuild(grid_size, grid_coord);
/**
 * figure of atom at z dimension if x is odd.
      z --->
      0.5    0   -0.5    -1    -1.5    -2   -2.5   -3   -3.5
      |____________|____________|____________|____________|
index:0            -1           -2          -3           -4
z lattice coord:
  z in [0, 1] -> index z = 0;
  z in [-1, 0] -> index z = -1;
  z in [-2, -1] -> index z = -2;
  z in [-3, -2] -> index z = -3;
 */
    // case 1: make z coord < 0.
    AtomElement src_atom;
    src_atom.x[0] = 3.5 * lattice_const; // odd x
    src_atom.x[1] = 0;
    src_atom.x[2] = 0 - 0.01;
    _type_atom_index coords[DIMENSION];
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[2], -1);

    // case 2:
    src_atom.x[2] = lattice_const * (-2) + 0.01;
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[2], -2);

    // case 3:
    src_atom.x[2] = lattice_const * (-2) - 0.01;
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[2], -3);

    // case 4:
    src_atom.x[2] = lattice_const * (-2.44);
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[2], -3);

    // case 5:
    src_atom.x[2] = lattice_const * (-2.66);
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[2], -3);
}

// test getNearLatSubBoxCoord when x[0] is minus.
TEST(ws_utils_getNearLatSubBoxCoord_minusX_test, _ws_utils_test) {
    const int grid_size[3] = {2, 2, 2};
    const int grid_coord[3] = {0, 0, 0};
    const int64_t space[3] = {25 * grid_size[0], 25 * grid_size[1], 25 * grid_size[2]};
    const double lattice_const = 2.85532;
    const double cutoff_radius_factor = 1.96125;
    comm::Domain *p_domain = comm::Domain::Builder()
            .setPhaseSpace(space)
            .setCutoffRadius(cutoff_radius_factor)
            .setLatticeConst(lattice_const)
            .localBuild(grid_size, grid_coord);

    // case 1: make x coord < 0.
    AtomElement src_atom;
    src_atom.x[0] = 0 - 0.01;
    _type_atom_index coords[DIMENSION];
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[0], 0 * 2);

    // case 2:
    src_atom.x[0] = lattice_const * (-2) - 0.01;
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[0], -2 * 2);

    // case 3:
    src_atom.x[0] = lattice_const * (-2 - 0.1); // point in (-1.75,-2.25) will be -2*2 = -4
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[0], -2 * 2);

    // case 4:
    src_atom.x[0] = lattice_const * (-2 + 0.1); // point in  (-1.75,-2.25) will be -2*2 = -4
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[0], -2 * 2);

    // case 5:
    src_atom.x[0] = lattice_const * (-2.25 - 0.1); // point in (-2.75,-2.25) will be -2.5*2 = -5
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[0], -2 * 2 - 1);

    // case 6:
    src_atom.x[0] = lattice_const * (-2.75 + 0.1); // point in (-2.75,-2.25) will be -2.5*2 = -5
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[0], -2 * 2 - 1);

    // case 7:
    src_atom.x[0] = lattice_const * (-2.75 - 0.1); // point in (-3,-2.75) will be -3*2 = -6
    ws::getNearLatSubBoxCoord(src_atom, p_domain, coords);
    EXPECT_EQ(coords[0], -3 * 2);
}
