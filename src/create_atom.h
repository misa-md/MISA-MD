//
// Created by baihe back to 2016-01-18.
//

#ifndef CRYSTAL_MD_CREATE_ATOM_H
#define CRYSTAL_MD_CREATE_ATOM_H

#include "atom.h"

#define boltz 8.617343e-5
#define mvv2e 1.0364269e-4

class create_atom {
public:
    create_atom(double tset);

    ~create_atom();

    void createphasespace(atom *_atom, double mass, int box_x, int box_y, int box_z);

private:
    double dof_compute(unsigned long natom);

    double t_set;
};

#endif // CRYSTAL_MD_CREATE_ATOM_H
