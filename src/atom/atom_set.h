//
// Created by genshen on 2018-12-31.
//

#ifndef CRYSTALMD_ATOM_SET_H
#define CRYSTALMD_ATOM_SET_H

#include <vector>

#include "domain/domain.h"
#include "atom/atom_element.h"
#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"

/**
 * the container of latticed atom list and inter atom list (all atoms in this sub-box).
 *
 */
class AtomSet {
public:
    AtomSet(Domain *domain, double latticeconst, double cutoffRadiusFactor);

    ~AtomSet();

    inline AtomList *getAtomList() {
        return atom_list;
    }

    inline AtomList &getAtomListRef() {
        return *atom_list;
    }

    inline InterAtomList *getInterList() {
        return inter_atom_list;
    }

    _type_atom_count getnlocalatom();

    /**
      * compute the index offset of neighbour atoms.
      */
    void calculateNeighbourIndices();

    /**
    * used in read creating mode.
    */
    void addAtom(unsigned long id, double rx, double ry, double rz, double vx, double vy, double vz);

protected:
    Domain *p_domain;

    int numberoflattice;

    double _cutoffRadius;
    int _cutlattice;
    double _latticeconst;

    std::vector<_type_atom_index> NeighbourOffsets; // 邻居粒子偏移量 // todo use offset in x,y,z dimension

    AtomList *atom_list;
    InterAtomList *inter_atom_list;
};


#endif //CRYSTALMD_ATOM_SET_H
