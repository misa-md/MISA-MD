//
// Created by genshen on 2019/11/13.
//

#ifndef CRYSTALMD_ARCH_IMP_H
#define CRYSTALMD_ARCH_IMP_H

#include <eam.h>
#include "atom/atom_element.h"
#include "arch_building_config.h"

#ifdef ACCELERATE_ENABLED

#define _ARCH_PREFIX_(arch_name, func_name)   arch_name ## _ ## func_name
#define ARCH_PREFIX(arch_name, func_name) _ARCH_PREFIX_(arch_name, func_name)

void ARCH_PREFIX(ARCH_NAME, env_init)();

void ARCH_PREFIX(ARCH_NAME, env_clean)();

void ARCH_PREFIX(ARCH_NAME, accelerate_init)(const comm::BccDomain *domain);

void ARCH_PREFIX(ARCH_NAME, pot_init)(eam *_pot);

void ARCH_PREFIX(ARCH_NAME, eam_rho_calc)(eam *pot, AtomElement *atoms, const double cutoff_radius);

void ARCH_PREFIX(ARCH_NAME, eam_df_calc)(eam *pot, AtomElement *atoms, const double cutoff_radius);

void ARCH_PREFIX(ARCH_NAME, eam_force_calc)(eam *pot, AtomElement *atoms, const double cutoff_radius);

#endif //ACCELERATE_ENABLED

#endif //CRYSTALMD_ARCH_IMP_H
