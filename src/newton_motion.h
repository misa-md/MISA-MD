//
// Created by genshen on 5/13/18.
//

#ifndef MISA_MD_NEWTON_H
#define MISA_MD_NEWTON_H

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

    void setTimestepLength(const double dt);

private:
    double _timestepLength;

    /**
     * index is consistent with the order in atom_type::atom_type in file atom_types.h
     */
    double dt_inv_m[atom_type::num_atom_types];

    void preComputeDtInv2m();
};


#endif //MISA_MD_NEWTON_H
