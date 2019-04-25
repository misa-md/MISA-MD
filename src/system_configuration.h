//
// Created by genshen on 2019-04-25.
//

#ifndef CRYSTALMD_SYSTEM_CONFIGURATION_H
#define CRYSTALMD_SYSTEM_CONFIGURATION_H


#include <array>

#include "types/pre_define.h"
#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"

/**
 * system configuration:
 * Api of statistical values of simulation system,
 * like temperature, system force, kinetic energy and potential energy.
 */
namespace configuration {
    /**
     * in development mode, we can check system temperature, energy or momentum.
     */
    std::array<_type_atom_force, DIMENSION> systemForce(AtomList *atom_list, InterAtomList *inter_atom_list);

};


#endif //CRYSTALMD_SYSTEM_CONFIGURATION_H
