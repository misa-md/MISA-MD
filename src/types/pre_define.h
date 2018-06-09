//
// Created by genshen on 5/6/18.
//

#ifndef CRYSTAL_MD_PRE_DEFINE_H
#define CRYSTAL_MD_PRE_DEFINE_H

#include "pre_config.h"

// constance
#define BOLTZ 8.617343e-5 // Boltzmann constant, 8.617343e-5 EV/K; also equals to 1.3806505e-23 J/K.
#define mvv2e 1.0364269e-4
#define ftm2v ( 1.0 / mvv2e)


// some value definition here.
#define DIMENSION 3

// domain
#define COORDINATE_ATOM_OUT_BOX (-100)

typedef int _type_lattice_size;
typedef _type_lattice_size _type_lattice_coord;

// atom
typedef unsigned long _type_atom_id;
typedef unsigned long _type_atom_index;
typedef unsigned long _type_atom_count;
typedef int _type_atom_type; // @deprecated
typedef short _type_inter_type; // inter atom type.
typedef unsigned short _type_atom_types; // the count of all types of atoms.

typedef double _type_atom_mass;
typedef double _type_atom_location;
typedef double _type_atom_velocity;
typedef double _type_atom_force;
typedef double _type_atom_rho;
typedef double _type_atom_df;

// simulation default default valuse
#include "../def_config_values.h"

#endif //CRYSTAL_MD_PRE_DEFINE_H
