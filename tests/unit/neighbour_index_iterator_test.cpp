//
// Created by genshen on 2019-04-14.
//

#include <gtest/gtest.h>
#include <types/pre_define.h>
#include <atom/neighbour_index.h>
#include <cmath>

class NeiIndexItlTests : public NeighbourIndex<AtomElement> {
public:
    explicit NeiIndexItlTests(AtomList &atom_list) : NeighbourIndex<AtomElement>(atom_list) {}

    FRIEND_TEST(nei_index_test_itl_plus, nei_index_test);
};

TEST(nei_index_test_itl_begin_end, nei_index_test) {
    const _type_lattice_size ext_size = 8;
    const _type_lattice_size box_size = 6;
    const _type_lattice_size ghost_size = ext_size - box_size;
    AtomList atom_list(ext_size * 2, ext_size, ext_size,
                       box_size * 2, box_size, box_size,
                       ghost_size * 2, ghost_size, ghost_size);
    NeiIndexItlTests nei_index(atom_list);
    // nei_index.make(ext_size - box_size, 1);

    // empty neighbour indexes vector
    NeiIndexItlTests::iterator begin = nei_index.begin(false, 0, 0, 0);
    NeiIndexItlTests::iterator end = nei_index.end(false, 0, 0, 0);
    if (begin != end) {
        GTEST_FAIL();
    }
}

TEST(nei_index_test_itl_plus, nei_index_test) {
    const _type_lattice_size ext_size = 8;
    const _type_lattice_size box_size = 6;
    const _type_lattice_size ghost_size = ext_size - box_size;
    AtomList atom_list(ext_size * 2, ext_size, ext_size,
                       box_size * 2, box_size, box_size,
                       ghost_size * 2, ghost_size, ghost_size);
    NeiIndexItlTests nei_index(atom_list);
    // nei_index.make(ext_size - box_size, 1);

    const int size_full_even = 2;
    const int size_full_odd = 4;
    const int size_half_even = 6;
    const int size_half_odd = 8;

    nei_index.nei_even_offsets.resize(size_full_even);
    nei_index.nei_odd_offsets.resize(size_full_odd);
    nei_index.nei_half_even_offsets.resize(size_half_even);
    nei_index.nei_half_odd_offsets.resize(size_half_odd);

    // full-even test
    {
        NeiIndexItlTests::iterator itl = nei_index.begin(false, 0, 0, 0);
        for (int i = 0; i < size_full_even; i++) {
            ++itl;
        }
        NeiIndexItlTests::iterator end = nei_index.end(false, 0, 0, 0);
        if (itl != end) {
            GTEST_FAIL();
        }
    }
    // full-odd test
    {
        NeiIndexItlTests::iterator itl = nei_index.begin(false, 3, 0, 0);
        for (int i = 0; i < size_full_odd; i++) {
            ++itl;
        }
        NeiIndexItlTests::iterator end = nei_index.end(false, 3, 0, 0);
        if (itl != end) {
            GTEST_FAIL();
        }
    }

    // half-even test
    {
        NeiIndexItlTests::iterator itl = nei_index.begin(true, 0, 0, 0);
        for (int i = 0; i < size_half_even; i++) {
            ++itl;
        }
        NeiIndexItlTests::iterator end = nei_index.end(true, 0, 0, 0);
        if (itl != end) {
            GTEST_FAIL();
        }
    }
    // half-odd test
    {
        NeiIndexItlTests::iterator itl = nei_index.begin(true, 3, 0, 0);
        for (int i = 0; i < size_half_odd; i++) {
            ++itl;
        }
        NeiIndexItlTests::iterator end = nei_index.end(true, 3, 0, 0);
        if (itl != end) {
            GTEST_FAIL();
        }
    }
}
