//
// Created by genshen on 2019-04-25.
//

#include "system_configuration.h"
#include "atom/atom_element.h"

std::array<_type_atom_force, DIMENSION> configuration::systemForce(
        AtomList *atom_list, InterAtomList *inter_atom_list) {
    _type_atom_force force_x = 0.0, force_y = 0.0, force_z = 0.0;
    atom_list->foreachSubBoxAtom([&force_x, &force_y, &force_z](AtomElement &_atom_ref) {
        if (_atom_ref.type != atom_type::INVALID) {
            force_x += _atom_ref.f[0];
            force_y += _atom_ref.f[1];
            force_z += _atom_ref.f[2];
        }
    });
    for (AtomElement &atom_ele:inter_atom_list->inter_list) {
        force_x += atom_ele.f[0];
        force_y += atom_ele.f[1];
        force_z += atom_ele.f[2];
    }
    return std::array<_type_atom_force, DIMENSION>{force_x, force_y, force_z};
}
