//
// Created by genshen on 2019-04-11.
//

#include "neighbour_index.h"


template<class T, class Ref, class Ptr>
NeiIterator<T, Ref, Ptr>::NeiIterator(): atom_list(nullptr) {}

template<class T, class Ref, class Ptr>
NeiIterator<T, Ref, Ptr>::NeiIterator(const NeiIterator::iterator &x)
        : atom_list(x.atom_list) {
    src_index = x.src_index;
    cur_index = x.cur_index;
    current_nei_index = x.current_nei_index;
    p_nei_index = x.p_nei_index;
}

template<class T, class Ref, class Ptr>
NeiIterator<T, Ref, Ptr>::NeiIterator(const link_type nei_index, const NeiOffset src_index,
                                      const AtomList *atom_list)
        :NeiIterator(nei_index, src_index, atom_list, 0) {}

template<class T, class Ref, class Ptr>
NeiIterator<T, Ref, Ptr>::NeiIterator(const link_type nei_index, const NeiOffset src_index,
                                      const AtomList *atom_list, NeiIterator::size_type offset)
        :p_nei_index(nei_index), src_index(src_index),
         cur_index(0), current_nei_index(offset), atom_list(atom_list) {
    // find first avail particle (some neighbours may be out of box) by
    // calculate particle index
    if (current_nei_index < p_nei_index->size()) {
        cur_index = src_index + (*p_nei_index)[current_nei_index];
    }
}

template<class T, class Ref, class Ptr>
NeiIterator<T, Ref, Ptr> &NeiIterator<T, Ref, Ptr>::operator++() {
    // find next avail neighbour particle.
    if (++current_nei_index < p_nei_index->size()) {
        // calculate particle index
        cur_index = src_index + (*p_nei_index)[current_nei_index];
        // note: with the ghost area, all neighbour particle can not be out of box.
    }
    return *this;
}
