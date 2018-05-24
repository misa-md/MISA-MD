//
// Created by genshen on 5/13/18.
//

#ifndef CRYSTAL_MD_NEWTON_H
#define CRYSTAL_MD_NEWTON_H

#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"

/**
 * Hat off to grate physicists.
 * resolving the equation of motion based on the Newton's Second Law.
 */
class NewtonMotion {
public:
    NewtonMotion(double timestepLength);

    ~NewtonMotion();

    void firststep(AtomList *atom_list, InterAtomList *inter_atom_list);

    void secondstep(AtomList *atom_list, InterAtomList *inter_atom_list);

    inline void setTimestepLength(double dt) {
        _timestepLength = dt;
    }

private:
    double _timestepLength;

    void computeFirst(double dtInv2m, AtomList *atom_list, InterAtomList *inter_atom_list);

    void computeSecond(double dtInv2m, AtomList *atom_list, InterAtomList *inter_atom_list);
};


#endif //CRYSTAL_MD_NEWTON_H
