//
// Created by genshen on 2021-06-16.
//

#include "atom_types.h"

std::vector<double> atom_type::mass_array;

_type_atom_types atom_type::num_atom_types;

void atom_type::setGlobalAtomMasses(const std::vector<double> masses) {
    num_atom_types = masses.size();
    mass_array.resize(masses.size());
    for (int i = 0; i < num_atom_types; i++) {
        mass_array[i] = masses[i];
    }
}
