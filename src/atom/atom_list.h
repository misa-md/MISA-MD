//
// Created by genshen on 5/8/18.
//

#ifndef MISA_MD_ATOM_LIST_H
#define MISA_MD_ATOM_LIST_H

#include <iterator>
#include <vector>
#include <functional>

#include <comm/domain/domain.h>
#include <comm/domain/bcc_domain.h>

#include "../types/pre_define.h"
#include "atom_element.h"
#include "lattice/lattice.h"


class AtomList {
public:
    friend class WorldBuilder;

    friend class atom;

    friend class AtomDump;

    // @see https://gist.github.com/jeetsukumaran/307264
    // @see https://stackoverflow.com/questions/3846296/how-to-overload-the-operator-in-two-different-ways-for-postfix-a-and-prefix
    // dont use it.
    class iterator {
    public:
        typedef iterator self_type;
        typedef std::forward_iterator_tag iterator_category;

        iterator(AtomElement *ptr);
//        iterator() {}
//
//        self_type &operator=(const self_type &other) { // todo copy.
//            ptr_ = other.ptr_;
//            return *this;
//        }

        self_type &operator++(); // PREFIX

        self_type operator++(int);  // POST

        AtomElement &operator*();

        AtomElement *operator->();

        bool operator==(const self_type &rhs);// todo non-const

        bool operator!=(const self_type &rhs);

    private:
        AtomElement *ptr_;

        void next();
    };

    class const_iterator {
    public:
        typedef const_iterator self_type;
        typedef std::forward_iterator_tag iterator_category;

        const_iterator(AtomElement *ptr);

        self_type &operator++();

        self_type operator++(int junk);

        const AtomElement &operator*();

        const AtomElement *operator->() { return ptr_; }

        bool operator==(const self_type &rhs);

        bool operator!=(const self_type &rhs);

    private:
        AtomElement *ptr_;
    };

    iterator begin();

    iterator end();

    const_iterator begin() const;

    const_iterator end() const;

    // todo override []
    // todo document.
    /**
     * initialize atom array with the size of atom(including ghost atoms).
     * @param size_x the atoms count at x dimension(including ghost atoms).
     * @param size_y the same as above at y dimension.
     * @param size_z the same as above at z dimension.
     * @param ghost_count_x the count of ghost atoms at left or right of x dimension. (in fact left and right has same count of ghost atoms)
     * @param ghost_count_y the same as above at y dimension.
     * @param ghost_count_z the same as above at z dimension.
     */
    AtomList(_type_atom_count size_x, _type_atom_count size_y, _type_atom_count size_z,
             _type_atom_count size_sub_box_x, _type_atom_count size_sub_box_y, _type_atom_count size_sub_box_z,
             _type_atom_count ghost_count_x,
             _type_atom_count ghost_count_y,
             _type_atom_count ghost_count_z);

    ~AtomList();

    /**
     * calculate index of 1d atom array by the lattice coordinate in 3d
     * @param x, y, z the lattice coordinate in 3d
     * @return the index in 1d atoms array of the atom specified by @param x,y,z coordinate
     */
    inline _type_atom_index getAtomIndex(_type_atom_index x, _type_atom_index y, _type_atom_index z) const {
        return (z * lattice._size_y + y) * lattice._size_x + x;
    }

    /**
     * get the reference of {@class AtomElement} by the lattice index in sub-box(not include ghost lattice).
     * @param index_x lattice index in sub-box at x dimension.
     * @param index_y lattice index in sub-box at y dimension.
     * @param index_z lattice index in sub-box at z dimension.
     * @return reference of {@class AtomElement} at corresponding position.
     */
    inline AtomElement &
    getAtomEleBySubBoxIndex(_type_atom_index index_x, _type_atom_index index_y, _type_atom_index index_z) {
        return _atoms[getAtomIndex(lattice.purge_ghost_count_x + index_x,
                                   lattice.purge_ghost_count_y + index_y,
                                   lattice.purge_ghost_count_z + index_z)];
    }

    /**
     * just by index of ghost lattice.
     * @param index_x
     * @param index_y
     * @param index_z
     * @return
     */
    inline AtomElement &
    getAtomEleByGhostIndex(_type_atom_index index_x, _type_atom_index index_y, _type_atom_index index_z) const {
        return _atoms[getAtomIndex(index_x, index_y, index_z)];
    }

    /**
     * get atom element by linear index.
     * index = (zIndex * p_domain->getGhostLatticeSize(1) + yIndex) * p_domain->getGhostLatticeSize(0) + xIndex;
     * @param index
     * @return
     */
    inline AtomElement &getAtomEleByLinearIndex(_type_atom_index index) const {
        return _atoms[index];
    }

    /**
     *
     */
    template<typename Callable>
    void foreachSubBoxAtom(Callable callback);

    /**
     * get the capacity of lattice atom in the box (including ghost lattice atoms).
     * or the lattice count (including ghost lattices).
     * @return
     */
    inline _type_lattice_size cap() const {
        return lattice._size;
    }

    /**
     * append/set an atom to inter. // todo unit test.
     * @param atom_id atom id.
     */
    void appendInter(_type_atom_id atom_id);

    void exchangeAtomFirst(comm::BccDomain *p_domain);

    void exchangeAtom(comm::BccDomain *p_domain);

    /**
     * return true if the there is atom in current box that is far away out of this box.
     * @param domain
     * @return
     */
    bool isBadList(comm::Domain domain);

public:
    const BccLattice lattice;

private:
    // 晶格点原子用数组存储其信息,including ghost atoms.
    AtomElement *_atoms; // atoms in 3d.

    // the array to record atoms that are out of box.
    std::vector<std::vector<_type_atom_id> > sendlist; // todo make it temp data
    std::vector<std::vector<_type_atom_id> > recvlist;

};


template<typename Callable>
void AtomList::foreachSubBoxAtom(Callable callback) {
    for (long z = lattice.purge_ghost_count_z; z < lattice._size_sub_box_z + lattice.purge_ghost_count_z; z++) {
        for (long y = lattice.purge_ghost_count_y; y < lattice._size_sub_box_y + lattice.purge_ghost_count_y; y++) {
            for (long x = lattice.purge_ghost_count_x; x < lattice._size_sub_box_x + lattice.purge_ghost_count_x; x++) {
                callback(getAtomEleByGhostIndex(x, y, z));
            }
        }
    }
}

#endif //MISA_MD_ATOM_LIST_H
