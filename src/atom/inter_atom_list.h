//
// Created by genshen on 2018-05-19.
//

#ifndef MISA_MD_INTER_ATOM_LIST_H
#define MISA_MD_INTER_ATOM_LIST_H

#include <vector>
#include <list>
#include <unordered_map>
#include <comm/domain/domain.h>
#include <comm/domain/bcc_domain.h>

#include "atom_list.h"
#include "atom_element.h"
#include "../types/pre_define.h"
#include "../types/atom_types.h"
#include "lattice/box.h"

typedef std::list<AtomElement> _type_inter_list;
typedef std::vector<std::vector<AtomElement *> > _type_inter_buf;
typedef std::pair<std::unordered_multimap<_type_atom_index, AtomElement *>::iterator,
        std::unordered_multimap<_type_atom_index, AtomElement *>::iterator> inter_map_range;
typedef std::unordered_multimap<_type_atom_index, AtomElement *>::iterator inter_map_range_itl;

/**
 * storing inter atoms
 */
class InterAtomList {
public:
    _type_inter_list inter_list;
    _type_inter_list inter_ghost_list;
    size_t nlocalinter; // 本地间隙原子数
    size_t nghostinter; // ghost间隙原子数

    std::unordered_multimap<_type_atom_index, AtomElement *> inter_map;

    InterAtomList();

    void appendInter(_type_atom_id atom_id);

    /**
     * add an inter atom to @var inter_list.
     * The atom data will be copied into the list.
     * @param atom reference of the atom.
     */
    void addInterAtom(AtomElement &atom);

    /**
     * insert an atom into ghost list, and return the atom pointer inserted.
     * @param ghost_atom ref of ghost atom.
     * @return the pointer inserted in atom list.
     */
    AtomElement *addGhostAtom(AtomElement &ghost_atom);

    inline size_t nLocalInter() {
        return nlocalinter;
    }

    void makeIndex(AtomList *atom_list, const comm::Domain *p_domain);

    _type_inter_list::iterator removeInter(_type_inter_list::iterator);

    void clearGhost();

};


#endif //MISA_MD_INTER_ATOM_LIST_H
