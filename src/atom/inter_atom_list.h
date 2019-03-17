//
// Created by genshen on 2018-05-19.
//

#ifndef CRYSTAL_MD_INTER_ATOM_LIST_H
#define CRYSTAL_MD_INTER_ATOM_LIST_H

#include <vector>
#include <list>

#include <domain/domain.h>

#include "atom_element.h"
#include "../pack/particledata.h"
#include "../pack/lat_particle_data.h"
#include "../types/pre_define.h"
#include "../types/atom_types.h"
#include "box.h"

typedef std::list<AtomElement> _type_inter_list;
typedef std::vector<std::vector<AtomElement *> > _type_inter_buf;

/**
 * storing inter atoms
 */
class InterAtomList {
public:
    _type_inter_list inter_list;
    _type_inter_list inter_ghost_list;
    size_t nlocalinter; // 本地间隙原子数
    size_t nghostinter; // ghost间隙原子数

    _type_inter_buf intersendlist;
    _type_inter_buf interrecvlist;

    InterAtomList();

    void appendInter(_type_atom_id atom_id);

    /**
     * add an inter atom to @var inter_list.
     * The atom data will be copied into the list.
     * @param atom reference of the atom.
     */
    void addInterAtom(AtomElement &atom);

    inline size_t nLocalInter() {
        return nlocalinter;
    }

    void exchangeInter(comm::Domain *p_domain);

    void borderInter(comm::Domain *p_domain);

    /**
     * pointer of element in atom_list (pointer of {@class AtomElement}).
     * // todo use avl tree.
     * // todo use pointer.
     */

};


#endif //CRYSTAL_MD_INTER_ATOM_LIST_H
