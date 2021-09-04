//
// Created by genshen on 2019-04-13.
//

#include <vector>
#include <array>
#include <gtest/gtest.h>
#include <types/pre_define.h>
#include <atom/neighbour_index.h>
#include <cmath>
#include <atom/atom_list.h>

class NeiIndexTests : public NeighbourIndex<_type_neighbour_index_ele> {
public:
    explicit NeiIndexTests(AtomList &atom_list)
            : NeighbourIndex<_type_neighbour_index_ele>(atom_list._atoms._data(), atom_list.lattice) {}

    FRIEND_TEST(isPositive_test, nei_index_test);

    FRIEND_TEST(nei_index_test_index_vector_len, nei_index_test);

    FRIEND_TEST(nei_index_test_index_1nn, nei_index_test);

    FRIEND_TEST(nei_index_test_index_2nn, nei_index_test);

    FRIEND_TEST(nei_index_test_index_case, nei_index_test);

    FRIEND_TEST(nei_index_test_compare_diff_methods, nei_index_test);

    // offset calculation in legacy CrystalMD
    std::array<std::vector<NeiOffset>, 4> legacyOffset(const _type_lattice_size cut_lattice,
                                                       const double cutoff_radius_factor) {
        std::vector<NeiOffset> _nei_even_offsets;
        std::vector<NeiOffset> _nei_odd_offsets;
        std::vector<NeiOffset> _nei_half_even_offsets;
        std::vector<NeiOffset> _nei_half_odd_offsets;
        const double cutoff_lat_factor = cutoff_radius_factor + 1.0 / 2 + 1.0 / 2;
        // if x index of a particle is even (the particle is lattice point,晶格点).
        for (_type_atom_index zIndex = -cut_lattice - 1;
             zIndex <= cut_lattice + 1; zIndex++) { // loop for (2*_cutlattice + 1) times.
            for (_type_atom_index yIndex = -cut_lattice - 1; yIndex <= cut_lattice + 1; yIndex++) {
                for (_type_atom_index xIndex = -cut_lattice - 1; xIndex <= cut_lattice + 1; xIndex++) {
                    // lattice neighbour points whose index is (2*Index,yIndex,zIndex).
                    {
                        const double r = xIndex * xIndex + yIndex * yIndex + zIndex * zIndex;
                        // r > 0 means neighbour index can not be itself.
                        if (r < cutoff_lat_factor * cutoff_lat_factor && r > 0) {
                            const _type_atom_index even_offset = lattice.IndexOf3DIndex(
                                    2 * xIndex, yIndex, zIndex);
                            _nei_even_offsets.push_back(even_offset);
                            if (isPositiveIndex(xIndex, yIndex, zIndex)) {
                                _nei_half_even_offsets.push_back(even_offset);
                            }
                        }
                    }
                    // bcc body center neighbour points whose index is (2*Index+1,yIndex,zIndex).
                    {
                        const double r = (xIndex + 0.5) * (xIndex + 0.5) +
                                         (yIndex + 0.5) * (yIndex + 0.5) +
                                         (zIndex + 0.5) * (zIndex + 0.5);
                        if (r < cutoff_lat_factor * cutoff_lat_factor && r > 0) { // in fact "r > 0" is not used.
                            const _type_atom_index even_offset = lattice.IndexOf3DIndex(
                                    2 * xIndex + 1, yIndex, zIndex);
                            _nei_even_offsets.push_back(even_offset);
                            if (isPositiveIndex(xIndex + 0.5, yIndex + 0.5, zIndex + 0.5)) {
                                _nei_half_even_offsets.push_back(even_offset);
                            }
                        }
                    }
                }
            }
        }
        // if x index of a particle is odd (the particle is BCC body center point,体心).
        for (_type_atom_index zIndex = -cut_lattice - 1;
             zIndex <= cut_lattice + 1; zIndex++) { // loop for (2*_cutlattice + 1) times.
            for (_type_atom_index yIndex = -cut_lattice - 1; yIndex <= cut_lattice + 1; yIndex++) {
                for (_type_atom_index xIndex = -cut_lattice - 1; xIndex <= cut_lattice + 1; xIndex++) {
                    // BCC body center neighbour points whose index is (2*Index,yIndex,zIndex).
                    {
                        const double r = xIndex * xIndex + yIndex * yIndex + zIndex * zIndex;
                        if (r < cutoff_lat_factor * cutoff_lat_factor && r > 0) {
                            const _type_atom_index odd_offset = lattice.IndexOf3DIndex(
                                    2 * xIndex, yIndex, zIndex);
                            _nei_odd_offsets.push_back(odd_offset);
                            if (isPositiveIndex(xIndex, yIndex, zIndex)) {
                                _nei_half_odd_offsets.push_back(odd_offset);
                            }
                        }
                    }
                    // lattice neighbour points whose index is (2*Index-1,yIndex,zIndex).
                    {
                        const double r = (xIndex - 0.5) * (xIndex - 0.5) +
                                         (yIndex - 0.5) * (yIndex - 0.5) +
                                         (zIndex - 0.5) * (zIndex - 0.5);
                        if (r < cutoff_lat_factor * cutoff_lat_factor && r > 0) {
                            const _type_atom_index odd_offset = lattice.IndexOf3DIndex(
                                    2 * xIndex - 1, yIndex, zIndex);
                            _nei_odd_offsets.push_back(odd_offset);
                            if (isPositiveIndex(xIndex - 0.5, yIndex - 0.5, zIndex - 0.6)) {
                                _nei_half_odd_offsets.push_back(odd_offset);
                            }
                        }
                    }
                }
            }
        }
        return std::array<std::vector<NeiOffset>, 4>{_nei_even_offsets, _nei_odd_offsets,
                                                     _nei_half_even_offsets, _nei_half_odd_offsets};
    }
};

TEST(isPositive_test, nei_index_test) {
    EXPECT_EQ(NeiIndexTests::isPositiveIndex(-3, -2, 0), false);
    EXPECT_EQ(NeiIndexTests::isPositiveIndex(-3, -2, -1), false);
    EXPECT_EQ(NeiIndexTests::isPositiveIndex(-3, 1, 0), true);
}

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

    EXPECT_EQ(nei_index.nei_odd_offsets.size(), 2 * nei_index.nei_half_odd_offsets.size());
    EXPECT_EQ(nei_index.nei_even_offsets.size(), 2 * nei_index.nei_half_even_offsets.size());
    EXPECT_EQ(nei_index.nei_odd_offsets.size(), nei_index.nei_even_offsets.size());
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

// fixme: wrong test case (this case will fail):
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
    NeighbourIndex<_type_neighbour_index_ele>::iterator nei_itl_end = nei_index.end(true, 24, 16, 3);
    for (NeighbourIndex<_type_neighbour_index_ele>::iterator nei_itl = nei_index.begin(true, 24, 16, 3);
         nei_itl != nei_itl_end; ++nei_itl) {
        _type_atom_index ind = atom_list._atoms.getAtomIndex(23, 17, 4);
        if (ind == nei_itl.cur_index) {
            exists = true;
        }
    }
    EXPECT_EQ(exists, true);
}

// compare offset index with another calculation method
TEST(nei_index_test_compare_diff_methods, nei_index_test) {
    const double cutoff_factor = 1.96125;
    const _type_lattice_size cutoff_lat = 2;
    const _type_lattice_size ghost_size = cutoff_lat;
    const _type_lattice_size box_size = 6;
    const _type_lattice_size ext_size = box_size + 2 * cutoff_lat;
    AtomList atom_list(ext_size * 2, ext_size, ext_size,
                       box_size * 2, box_size, box_size,
                       ghost_size * 2, ghost_size, ghost_size);
    NeiIndexTests nei_index(atom_list);
    nei_index.make(cutoff_lat, cutoff_factor);
    auto offsets_array = nei_index.legacyOffset(cutoff_lat, cutoff_factor);

    auto legacy_nei_even_offsets = offsets_array[0];
    auto legacy_nei_odd_offsets = offsets_array[1];
    auto legacy_nei_half_even_offsets = offsets_array[2];
    auto legacy_nei_half_odd_offsets = offsets_array[3];

    // compare size
    EXPECT_EQ(legacy_nei_even_offsets.size(), nei_index.nei_even_offsets.size());
    EXPECT_EQ(legacy_nei_odd_offsets.size(), nei_index.nei_odd_offsets.size());
    EXPECT_EQ(legacy_nei_half_even_offsets.size(), nei_index.nei_half_even_offsets.size());
    EXPECT_EQ(legacy_nei_half_odd_offsets.size(), nei_index.nei_half_odd_offsets.size());

    // sort and compare elements.
    std::sort(legacy_nei_even_offsets.begin(), legacy_nei_even_offsets.end());
    std::sort(legacy_nei_odd_offsets.begin(), legacy_nei_odd_offsets.end());
    std::sort(legacy_nei_half_even_offsets.begin(), legacy_nei_half_even_offsets.end());
    std::sort(legacy_nei_half_odd_offsets.begin(), legacy_nei_half_odd_offsets.end());

    std::sort(nei_index.nei_even_offsets.begin(), nei_index.nei_even_offsets.end());
    std::sort(nei_index.nei_odd_offsets.begin(), nei_index.nei_odd_offsets.end());
    std::sort(nei_index.nei_half_even_offsets.begin(), nei_index.nei_half_even_offsets.end());
    std::sort(nei_index.nei_half_odd_offsets.begin(), nei_index.nei_half_odd_offsets.end());

    size_t i = 0;
    for (auto ele:legacy_nei_even_offsets) {
        EXPECT_EQ(ele, nei_index.nei_even_offsets[i]);
        i++;
    }
    i = 0;
    for (auto ele:legacy_nei_odd_offsets) {
        EXPECT_EQ(ele, nei_index.nei_odd_offsets[i]);
        i++;
    }
    i = 0;
    for (auto ele:legacy_nei_half_even_offsets) {
        EXPECT_EQ(ele, nei_index.nei_half_even_offsets[i]);
        i++;
    }
    i = 0;
    for (auto ele:legacy_nei_half_odd_offsets) {
        EXPECT_EQ(ele, nei_index.nei_half_odd_offsets[i]);
        i++;
    }
}
