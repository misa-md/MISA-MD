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

double configuration::kineticEnergy(AtomList *atom_list, InterAtomList *inter_atom_list,
                                    ReturnMod mode, const kiwi::RID root) {
    return kineticEnergy(atom_list, mode, root) + kineticEnergy(inter_atom_list, mode, root);
}

double configuration::kineticEnergy(AtomList *atom_list, ReturnMod mode, const kiwi::RID root) {
    double energy = 0;
    atom_list->foreachSubBoxAtom([&energy](AtomElement &_atom_ref) {
        if (_atom_ref.type != atom_type::INVALID) {
            // energy += v^2*m
            energy += (_atom_ref.v[0] * _atom_ref.v[0] +
                       _atom_ref.v[1] * _atom_ref.v[1] +
                       _atom_ref.v[2] * _atom_ref.v[2]) *
                      atom_type::getAtomMass(_atom_ref.type);
        }
    });
    // sum energy from all processes.
    double e_global = 0;
    switch (mode) {
        case Local:
            return energy / 2;
        case Root:
            MPI_Reduce(&energy, &e_global, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
            break;
        case All:
            MPI_Allreduce(&energy, &e_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            break;
    }
    return e_global / 2;
}

double configuration::kineticEnergy(InterAtomList *inter_atom_list, ReturnMod mode, const kiwi::RID root) {
    double energy = 0;
    for (_type_inter_list::iterator itl = inter_atom_list->inter_list.begin();
         itl != inter_atom_list->inter_list.end(); ++itl) {
        AtomElement &_atom_ref = *itl;
        energy += (_atom_ref.v[0] * _atom_ref.v[0] +
                   _atom_ref.v[1] * _atom_ref.v[1] +
                   _atom_ref.v[2] * _atom_ref.v[2]) *
                  atom_type::getAtomMass(_atom_ref.type);
    }
    // sum energy from all processes.
    double e_global = 0;
    switch (mode) {
        case Local:
            return energy / 2;
        case Root:
            MPI_Reduce(&energy, &e_global, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
            break;
        case All:
            MPI_Allreduce(&energy, &e_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            break;
    }
    return e_global / 2;
}

double configuration::temperature(const double ke, const _type_lattice_size n) {
//    dof -= 3; // fixme, why?
    // ke = 3nkT/2
    const double to_T = mvv2e / ((3 * n - 3) * BOLTZ); // 2 * ke * mvv2e / ((3 * n - 3) * BOLTZ); /
    return 2 * ke * to_T; // todo better times order for precision.
}
