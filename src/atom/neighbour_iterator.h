//
// Created by genshen on 2019-04-11.
//

#ifndef CRYSTAL_MD_NEIGHBOUR_ITERATOR_H
#define CRYSTAL_MD_NEIGHBOUR_ITERATOR_H

#include <iterator>
#include <cstddef>
#include <types/pre_define.h>
#include "atom_list.h"

template<typename T>
class NeighbourIndex;

template<class T, class Ref, class Ptr>
class NeiIterator {
public:
    typedef NeiIterator<T, T &, T *> iterator;
    typedef NeiIterator<T, const T &, const T *> const_iterator;
    typedef NeiIterator<T, Ref, Ptr> self;

    typedef std::forward_iterator_tag iterator_category;

    typedef T value_type;
    typedef Ptr pointer;
    typedef Ref reference;
    typedef std::vector<NeiOffset> *link_type;
//    typedef typename Bucket<T>::bucket_iterator bucket_iterator_type; // iterator of list in Bucket
    typedef size_t size_type;

    explicit NeiIterator();

    /**
     * iterator copy.
     */
    NeiIterator(const iterator &x);

    /**
     * \brief initialize iterator using neighbour index and source atom indexes.
     * set the iterator to the first element in neighbour indexes list.
     * \param nei_index neighbour indexes list.
     * \param src_index source index
     * \param atom_list atom list, all atoms here.
     */
    NeiIterator(const link_type nei_index, const NeiOffset src_index, const AtomList *atom_list);

    /**
     * \brief initialize iterator using neighbour indexes list, source index and index in neighbour indexes list.
     * \param nei_index neighbour indexes list.
     * \param src_index source index.
     * \param atom_list atom list, all atoms here.
     * \param offset offset for neighbour index list.
     */
    NeiIterator(const link_type nei_index, const NeiOffset src_index,
                const AtomList *atom_list, size_type offset);

    /**
     * iterator equal is specified by neighbour lattice indexes vector, atom list, current index offset and source lattice.
     * @param x another iterator
     * @return true: equal, false: not equal.
     */
    bool operator==(const self &x) const {
        return current_nei_index == x.current_nei_index &&
               p_nei_index == x.p_nei_index &&
               atom_list == x.atom_list &&
               src_index == x.src_index;
    }

    bool operator!=(const self &x) const {
        return !(*this == x);
    }

    reference operator*() {
        return atom_list->getAtomEleByLinearIndex(cur_index);
    }

    pointer operator->() { return &(operator*()); }

    // return reference of next available node.
    self &operator++();

public:
    // current index of particle, cur_index = src_index + p_nei_index[current_nei_index].
    _type_atom_index cur_index;

protected:
    const link_type p_nei_index; // pointer of neighbour index vector
    // source particle index for iterator, all offset of neighbour particles is based on this source particles index.
    const _type_atom_index src_index;
    // current index for neighbour index.
    size_type current_nei_index;
    const AtomList *atom_list;
};

#include "neighbour_iterator.inl"


#endif //CRYSTAL_MD_NEIGHBOUR_ITERATOR_H
