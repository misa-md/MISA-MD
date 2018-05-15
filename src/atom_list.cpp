//
// Created by genshen on 5/8/18.
//

#include "atom_list.h"

AtomList::AtomList(_type_atom_count size_x, _type_atom_count size_y, _type_atom_count size_z,
                   _type_atom_count ghost_count_x, _type_atom_count ghost_count_y, _type_atom_count ghost_count_z) :
        _size(size_x * size_y * size_z),
        _size_x(size_x), _size_y(size_y), _size_z(size_z),
        purge_ghost_count_x(ghost_count_x),
        purge_ghost_count_y(ghost_count_y),
        purge_ghost_count_z(ghost_count_z) {

    _atoms = new AtomElement **[size_z];
    for (_type_atom_count z = 0; z < size_z; z++) {
        _atoms[z] = new AtomElement *[size_y];
        for (_type_atom_count y = 0; y < size_y; y++) {
            _atoms[z][y] = new AtomElement[size_x];
        }
    }
}

AtomList::~AtomList() {
    for (_type_atom_count z = 0; z < _size_z; z++) {
        for (_type_atom_count y = 0; y < _size_y; y++) {
            delete[] _atoms[z][y];
        }
        delete[] _atoms[z];
    }
    delete[] _atoms;
}

void AtomList::appendInter(_type_atom_id atom_id) {

}

// iterator
AtomList::iterator::iterator(AtomElement *ptr) : ptr_(ptr) {}

AtomList::iterator::self_type &AtomList::iterator::operator++() {
    ptr_++;
    return *this;
}

AtomList::iterator::self_type AtomList::iterator::operator++(int) {
    self_type i = *this;
    ptr_++;
    return i;
}

AtomElement &AtomList::iterator::operator*() {
    return *ptr_;
}

AtomElement *AtomList::iterator::operator->() {
    return ptr_;
}

bool AtomList::iterator::operator==(const AtomList::iterator::self_type &rhs) {
    return ptr_ == rhs.ptr_;
}

bool AtomList::iterator::operator!=(const AtomList::iterator::self_type &rhs) {
    return ptr_ != rhs.ptr_;
}

void AtomList::iterator::next() {

}

// const iterator
AtomList::const_iterator::const_iterator(AtomElement *ptr) : ptr_(ptr) {}

AtomList::const_iterator::self_type &AtomList::const_iterator::operator++() {
    ptr_++;
    return *this;
}

AtomList::const_iterator::self_type AtomList::const_iterator::operator++(int junk) {
    self_type i = *this;
    ptr_++;
    return i;
}

const AtomElement &AtomList::const_iterator::operator*() {
    return *ptr_;
}

bool AtomList::const_iterator::operator==(const AtomList::const_iterator::self_type &rhs) {
    return ptr_ == rhs.ptr_;
}

bool AtomList::const_iterator::operator!=(const AtomList::const_iterator::self_type &rhs) {
    return ptr_ != rhs.ptr_;
}
