//
// Created by genshen on 5/6/18.
//

#ifndef CRYSTAL_MD_PRE_DEFINE_H
#define CRYSTAL_MD_PRE_DEFINE_H

#include "../md_building_config.h"

// constance
#define BOLTZ 8.617343e-5 // Boltzmann constant, 8.617343e-5 EV/K; also equals to 1.3806505e-23 J/K.

// x ev= x * 1.602176634×10−19 J = mvv/2
// => v' = sqrt(x * 2*1.602176634×10−19/m)  note: the unit of m is kg, unit of v' is meter/s. we have m = M*1.66053886e-27 kg-1
// => v' = sqrt(x * 2*C/M), where C = 1.602176634e-19/1.66053886e-27,
// let v in 100m/s (also A/ps), v = v'/100 = sqrt(x * 2 * C / M/10000).
// let C/10000 = 1/mvv2e, we can get: mvv2e = 10000/C = 0.000103642684
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
typedef long _type_atom_index; // index can be minus
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


#endif //CRYSTAL_MD_PRE_DEFINE_H
