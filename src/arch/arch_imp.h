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

void ARCH_PREFIX(ARCH_NAME, accelerate_init)(const int lolocalx, const int lolocaly, const int lolocalz,
                                             const int nlocalx, const int nlocaly, const int nlocalz,
                                             const int loghostx, const int loghosty, const int loghostz,
                                             const int nghostx, const int nghosty, const int nghostz);

void ARCH_PREFIX(ARCH_NAME, eam_rho_calc)(int *rho_n, AtomElement ***atoms, double *cutoffRadius,
                                          double *rhoInvDx, double *rhoSplineValues);

void ARCH_PREFIX(ARCH_NAME, eam_df_calc)(int *df_n, AtomElement ***atoms, double *cutoffRadius,
                                         double *dfSplineInvDx, double *dfSplineValues);

void ARCH_PREFIX(ARCH_NAME, eam_force_calc)(int *phi_n, AtomElement ***atoms,
                                            double *cutoffRadius, double *phiSplineInvDx,
                                            double *phiSplineValues, double *rhoSplineValues);

#endif //ACCELERATE_ENABLED

#endif //CRYSTALMD_ARCH_IMP_H
