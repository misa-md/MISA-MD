//
// Created by genshen on 5/7/18.
//

#ifndef CRYSTAL_MD_ATOM_ELEMENT_H
#define CRYSTAL_MD_ATOM_ELEMENT_H


#include "pre_config.h"

/**
 * This class describes the attributes of one atom, such as location, force, velocity.etc.
 */

typedef unsigned long _type_atom_id;
typedef int _type_atom_type;
typedef double _type_atom_location;
typedef double _type_atom_velocity;
typedef double _type_atom_force;
typedef double _type_atom_rho;
typedef double _type_atom_df;

class AtomElement {
public:
    _type_atom_id id; // atom id.
    _type_atom_type type; // atom type

    _type_atom_location x[DIMENSION]; // atom position.
    _type_atom_velocity v[DIMENSION]; // atom velocity.
    _type_atom_force f[DIMENSION];  // atom force.

    _type_atom_rho rho;
    _type_atom_rho df;
};


#endif //CRYSTAL_MD_ATOM_ELEMENT_H
