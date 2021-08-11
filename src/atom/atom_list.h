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
#include "atom_prop_list.hpp"
#include "lattice/lattice.h"


class AtomList {
public:
    friend class WorldBuilder;

    friend class atom;

    friend class AtomDump;

    // @see https://gist.github.com/jeetsukumaran/307264
    // @see https://stackoverflow.com/questions/3846296/how-to-overload-the-operator-in-two-different-ways-for-postfix-a-and-prefix
    // don't use it.
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
     *
     */
    template<typename Callable>
    void foreachSubBoxAtom(Callable callback);

    template<typename Callback>
    void foreachSubBoxAtoms(Callback callback);

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

    /**
     * return true if the there is atom in current box that is far away out of this box.
     * @param domain
     * @return
     */
    bool isBadList(comm::Domain domain);

public:
    const BccLattice lattice;

public:
    // 晶格点原子用数组存储其信息,including ghost atoms.
    AtomPropList<AtomElement> _atoms; // atoms in 3d.
};


template<typename Callable>
void AtomList::foreachSubBoxAtom(Callable callback) {
    for (long z = lattice.purge_ghost_count_z; z < lattice._size_sub_box_z + lattice.purge_ghost_count_z; z++) {
        for (long y = lattice.purge_ghost_count_y; y < lattice._size_sub_box_y + lattice.purge_ghost_count_y; y++) {
            for (long x = lattice.purge_ghost_count_x; x < lattice._size_sub_box_x + lattice.purge_ghost_count_x; x++) {
                callback(_atoms.getAtomEleByGhostIndex(x, y, z));
            }
        }
    }
}

#endif //MISA_MD_ATOM_LIST_H
