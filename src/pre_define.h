//
// Created by genshen on 5/6/18.
//

#ifndef CRYSTAL_MD_PRE_DEFINE_H
#define CRYSTAL_MD_PRE_DEFINE_H

#include "pre_config.h"

// constance
#define BOLTZ 8.617343e-5 // Boltzmann constant, 8.617343e-5 EV/K; also equals to 1.3806505e-23 J/K.
#define mvv2e 1.0364269e-4 // todo move to predefine.h


// some value definition here.
#define DIMENSION 3

// simulation
#define DEFAULT_TIME_STEP_LENGTH 0.001

// domain
#define COORDINATE_ATOM_OUT_BOX (-100)

typedef int _type_lattice_size;
typedef _type_lattice_size _type_lattice_coord;

// atom
typedef unsigned long _type_atom_id;

#endif //CRYSTAL_MD_PRE_DEFINE_H
