//
// Created by genshen on 2019-01-02.
//

#include <gtest/gtest.h>
#include <domain/domain.h>
#include <atom/inter_atom_list.h>
#include <atom/ws_utils.h>

Domain *buildLocalDomain(const int64_t local_space[DIMENSION],
                         const int grid_size[DIMENSION],
                         const int grid_coord[DIMENSION],
                         const double lattice_const = 0.86,
                         const double cutoff_radius_factor = 1.1421) {

    const int64_t space[3] = {local_space[0] * grid_size[0],
                              local_space[1] * grid_size[1],
                              local_space[2] * grid_size[2]};
    return Domain::Builder()
            .setPhaseSpace(space)
            .setCutoffRadius(cutoff_radius_factor)
            .setLatticeConst(lattice_const)
            .localBuild(grid_size, grid_coord);
}

TEST(inter_atom_test_isInBox_litter, inter_atom_Test) {
    const double lattice_const = 0.86;
    const double cutoff_factor = 1.1421;
    const int grid_size[3] = {2, 2, 2};
    const int grid_coord[3] = {0, 0, 0};
    const int64_t local_space[3] = {50, 60, 72};
    Domain *p_domain = buildLocalDomain(local_space, grid_size, grid_coord, lattice_const, cutoff_factor);

    AtomElement atom;
    // test in box
    atom.x[0] = 10 * lattice_const;
    atom.x[1] = 10 * lattice_const;
    atom.x[2] = 10 * lattice_const;
    auto flag1 = ws::isOutBox(atom, p_domain);
    EXPECT_EQ(flag1, box::IN_BOX);

    // test out box: x direction with litter end.
    atom.x[0] = -10 * lattice_const;
    atom.x[1] = 10 * lattice_const;
    atom.x[2] = 10 * lattice_const;
    auto flag2 = ws::isOutBox(atom, p_domain);
    EXPECT_EQ(flag2, box::OUT_BOX_X_LITTER);

    atom.x[0] = -0.01; // not out of box (the near atom is still on box).
    atom.x[1] = 10 * lattice_const;
    atom.x[2] = 10 * lattice_const;
    auto flag3 = ws::isOutBox(atom, p_domain);
    EXPECT_EQ(flag3, box::IN_BOX);

    atom.x[0] = -lattice_const / 2 - 0.01; // not out of box (the near atom is still on box).
    atom.x[1] = -lattice_const / 2 - 0.01;
    atom.x[2] = -lattice_const / 2 - 0.01;
    auto flag4 = ws::isOutBox(atom, p_domain);
    EXPECT_EQ(flag4, box::OUT_BOX_X_LITTER | box::OUT_BOX_Y_LITTER | box::OUT_BOX_Z_LITTER);

    delete p_domain;
}
