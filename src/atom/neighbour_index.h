//
// Created by genshen on 2019-04-11.
//

#ifndef CRYSTAL_MD_NEIGHBOUR_INDEX_H
#define CRYSTAL_MD_NEIGHBOUR_INDEX_H

#include <domain/region.hpp>
#include "neighbour_iterator.h"
#include "atom_list.h"

/**
 * NeighbourIndex stores the relative index of atoms in BCC box.
 * @tparam T type of particles
 */
template<class T>
class NeighbourIndex {
    friend class iterator;

public:
    typedef NeiIterator<T, T &, T *> iterator;

    /**
     * Initialize neighbour particles index.
     * @param atom_list ref of atom list.
     */
    explicit NeighbourIndex(AtomList &atom_list);

    /**
     * create neighbour index list.
     * @param cut_lattice the cutoff radius of neighbour lattices in lattice size.
     * @param cutoff_radius_factor the cutoff radius, it muse be less then @param cut_lattice.
     * @note the size_x is normal lattice size, which is not doubled due to the data structure..
     */
    void make(const _type_lattice_size cut_lattice, const double cutoff_radius_factor);

    /**
     * begin of iterator
     * @param half_itl whether to use half iterator(iterator nei_half_even_offsets or nei_half_odd_offsets).
     * @param x, y, z the coord of source lattices at each dimension.
     * @note x is twice times as the normal lattice size at x dimension.
     */
    iterator begin(const bool half_itl, const _type_atom_index x,
                   const _type_atom_index y, const _type_atom_index z);

    /**
     * end of iterator
     * @param half_itl whether to use half iterator.
     * @param x, y, z the coord of source lattices at each dimension.
     */
    iterator end(const bool half_itl, const _type_atom_index x,
                 const _type_atom_index y, const _type_atom_index z);

protected:
    const AtomList &atom_list;

    // particles offset indexes of neighbours lattices.
    std::vector<NeiOffset> nei_even_offsets;
    std::vector<NeiOffset> nei_odd_offsets;
    std::vector<NeiOffset> nei_half_even_offsets;
    std::vector<NeiOffset> nei_half_odd_offsets;

    /**
     * If this index is a positive index, true will be returned.
     * @param x,y,z the index in x,y,z dimension.
     * @return true for positive index, false for otherwise.
     */
    static bool isPositiveIndex(const double x, const double y, const double z);
};


#include "neighbour_index.inl"


#endif //CRYSTAL_MD_NEIGHBOUR_INDEX_H
