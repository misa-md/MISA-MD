//
// Created by genshen on 2019-04-11.
//

#include "neighbour_index.h"


template<class T, class Ref, class Ptr, class DataTp>
NeiIterator<T, Ref, Ptr, DataTp>::NeiIterator(): atom_list(nullptr) {}

template<class T, class Ref, class Ptr, class DataTp>
NeiIterator<T, Ref, Ptr, DataTp>::NeiIterator(const NeiIterator::iterator &x)
        :cur_index(x.cur_index),
         p_nei_index(x.p_nei_index),
         src_index(x.src_index),
         current_nei_index(x.current_nei_index),
         atom_list(x.atom_list) {
}

template<class T, class Ref, class Ptr, class DataTp>
NeiIterator<T, Ref, Ptr, DataTp>::NeiIterator(const link_type nei_index, const NeiOffset src_index,
                                              DataTp atom_list)
        :NeiIterator<T, Ref, Ptr, DataTp>(nei_index, src_index, atom_list, 0) {}

template<class T, class Ref, class Ptr, class DataTp>
NeiIterator<T, Ref, Ptr, DataTp>::NeiIterator(const link_type nei_index, const NeiOffset src_index,
                                              DataTp atom_list, NeiIterator::size_type offset)
        : cur_index(0), current_nei_index(offset),
          p_nei_index(nei_index), src_index(src_index), atom_list(atom_list) {
    // find first avail particle (some neighbours may be out of box) by
    // calculate particle index
    if (current_nei_index < p_nei_index->size()) {
        cur_index = src_index + (*p_nei_index)[current_nei_index];
    }
}

template<class T, class Ref, class Ptr, class DataTp>
NeiIterator<T, Ref, Ptr, DataTp> &NeiIterator<T, Ref, Ptr, DataTp>::operator++() {
    // find next avail neighbour particle.
    if (++current_nei_index < p_nei_index->size()) {
        // calculate particle index
        cur_index = src_index + (*p_nei_index)[current_nei_index];
        // note: with the ghost area, all neighbour particle can not be out of box.
    }
    return *this;
}
