//
// Created by genshen on 2019-04-13.
//

#include <gtest/gtest.h>
#include <types/pre_define.h>
#include <atom/neighbour_index.h>
#include <cmath>

class NeiIndexTests : public NeighbourIndex<AtomElement> {
public:
    explicit NeiIndexTests(AtomList &atom_list) : NeighbourIndex<AtomElement>(atom_list) {}

    FRIEND_TEST(nei_index_test_index_vector_len, nei_index_test);

    FRIEND_TEST(nei_index_test_index_1nn, nei_index_test);

    FRIEND_TEST(nei_index_test_index_2nn, nei_index_test);

    FRIEND_TEST(nei_index_test_index_case, nei_index_test);
};

// test index array length
TEST(nei_index_test_index_vector_len, nei_index_test) {
    const _type_lattice_size ext_size = 8;
    const _type_lattice_size box_size = 6;
    const _type_lattice_size ghost_size = ext_size - box_size;
    AtomList atom_list(ext_size * 2, ext_size, ext_size,
                       box_size * 2, box_size, box_size,
                       ghost_size * 2, ghost_size, ghost_size);
    NeiIndexTests nei_index(atom_list);
    nei_index.make(ghost_size, ghost_size - 0.1);
    EXPECT_EQ(nei_index.nei_odd_offsets.size(), 2 * nei_index.nei_half_even_offsets.size());
    EXPECT_EQ(nei_index.nei_even_offsets.size(), 2 * nei_index.nei_half_even_offsets.size());
}

TEST(nei_index_test_index_1nn, nei_index_test) {
    const _type_lattice_size ext_size = 8;
    const _type_lattice_size box_size = 6;
    const _type_lattice_size ghost_size = ext_size - box_size;
    AtomList atom_list(ext_size * 2, ext_size, ext_size,
                       box_size * 2, box_size, box_size,
                       ghost_size * 2, ghost_size, ghost_size);
    NeiIndexTests nei_index(atom_list);
    nei_index.make(2, 0.51 * sqrt(3)); // search 2 lattice size to find 1nn neighbours.
    // 8 1nn neighbours
    EXPECT_EQ(nei_index.nei_odd_offsets.size(), 8);
    EXPECT_EQ(nei_index.nei_even_offsets.size(), 8);
    // todo test more detail
}

TEST(nei_index_test_index_2nn, nei_index_test) {
    const _type_lattice_size ext_size = 8;
    const _type_lattice_size box_size = 6;
    const _type_lattice_size ghost_size = ext_size - box_size;
    AtomList atom_list(ext_size * 2, ext_size, ext_size,
                       box_size * 2, box_size, box_size,
                       ghost_size * 2, ghost_size, ghost_size);
    NeiIndexTests nei_index(atom_list);
    nei_index.make(2, 1.1 * sqrt(3)); // search 2 lattice size to find 1nn and 2nn neighbours.
    // 8 1nn neighbours
    EXPECT_EQ(nei_index.nei_odd_offsets.size(), 8 + 6);
    EXPECT_EQ(nei_index.nei_even_offsets.size(), 8 + 6);
    // todo test more detail
}

TEST(nei_index_test_index_case, nei_index_test) {
    const _type_lattice_size ext_size = 9;
    const _type_lattice_size box_size = 6;
    const _type_lattice_size ghost_size = ext_size - box_size;
    AtomList atom_list(ext_size * 2, ext_size, ext_size,
                       box_size * 2, box_size, box_size,
                       ghost_size * 3, ghost_size, ghost_size);
    NeiIndexTests nei_index(atom_list);
    nei_index.make(3, 1.96125);

//    24	16	3
// its neighbour (23,17,4) should be included.
    bool exists = false;
    NeighbourIndex<AtomElement>::iterator nei_itl_end = nei_index.end(true, 24, 16, 3);
    for (NeighbourIndex<AtomElement>::iterator nei_itl = nei_index.begin(true, 24, 16, 3);
         nei_itl != nei_itl_end; ++nei_itl) {
        if (nei_itl.cur_index_x == 23 && nei_itl.cur_index_y == 17 && nei_itl.cur_index_z == 4) {
            exists = true;
        }
    }
    EXPECT_EQ(exists, true);
}
