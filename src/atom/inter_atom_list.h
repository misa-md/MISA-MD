//
// Created by genshen on 2018-05-19.
//

#ifndef CRYSTAL_MD_INTER_ATOM_LIST_H
#define CRYSTAL_MD_INTER_ATOM_LIST_H

#include <vector>
#include <list>
#include <pack/lat_particle_data.h>
#include <pack/particledata.h>
#include "atom_element.h"
#include "../types/pre_define.h"
#include "../types/atom_types.h"

typedef std::list<AtomElement> _type_inter_list;

/**
 * storing inter atoms
 */
class InterAtomList {
public:

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

    void pack_intersend(std::vector<unsigned long> interbuf, particledata *buf);

    void unpack_interrecv(int d, int n,
                          double lower[DIMENSION], // p_domain->getMeasuredSubBoxLowerBounding(d)
                          double upper[DIMENSION], // p_domain->getMeasuredSubBoxUpperBounding(d)
                          particledata *buf);

    void pack_bordersend(int dimension, int n, std::vector<int> &sendlist, LatParticleData *buf, double shift);

    void unpack_borderrecv(int n, double lower[DIMENSION], // p_domain->getMeasuredGhostLowerBounding(d)
                           double upper[DIMENSION], // p_domain->getMeasuredGhostUpperBounding(d)
                           LatParticleData *buf, std::vector<int> &recvlist);

    _type_inter_list inter_list;
    _type_inter_list inter_ghost_list;
    size_t nlocalinter; // 本地间隙原子数
    size_t nghostinter; // ghost间隙原子数

    /**
     * pointer of element in atom_list (pointer of {@class AtomElement}).
     * // todo use avl tree.
     * // todo use pointer.
     */
};


#endif //CRYSTAL_MD_INTER_ATOM_LIST_H
