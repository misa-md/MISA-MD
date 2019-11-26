//
// Created by genshen on 2018-3-5.
//

#ifndef CRYSTAL_MD_HARDWARE_ACCELERATE_H
#define CRYSTAL_MD_HARDWARE_ACCELERATE_H

#include <eam.h>
#include "atom/atom_element.h"
#include "arch_building_config.h"

#include "arch_imp.h"

// check whether it has accelerate hardware to be used, for example GPU, MIC(Xeon Phi), or sunway slave cores.
inline bool isAccelerateSupport() {
#ifdef ACCELERATE_ENABLED
    return true; // sunway and other hardware.
#else
    return false;
#endif
}

// initial for hardware accelerate.
// about const &, see: https://stackoverflow.com/questions/9637856/why-is-const-int-faster-than-const-int/9637951#9637951
inline void accelerateInit(const int lolocalx, const int lolocaly, const int lolocalz,
                           const int nlocalx, const int nlocaly, const int nlocalz,
                           const int loghostx, const int loghosty, const int loghostz,
                           const int nghostx, const int nghosty, const int nghostz) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, accelerate_init)(lolocalx, lolocaly, lolocalz, nlocalx, nlocaly, nlocalz,
                                            loghostx, loghosty, loghostz, nghostx, nghosty, nghostz);
#endif
}

// it runs after atom and boxes creation, but before simulation running.
inline void beforeAccelerateRun(eam *_pot) {
#ifdef ACCELERATE_ENABLED
    // ARCH_PREFIX(ARCH_NAME, pot_init)(_pot->electron_density->spline, _pot->f->spline, _pot->phi->spline);
#endif
}

// accelerate for calculating electron_density in computing eam potential.
inline void accelerateEamRhoCalc(int *rho_n, AtomElement ***atoms, double *cutoffRadius,
                                 double *rhoInvDx, double *rhoSplineValues) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, eam_rho_calc)(rho_n, atoms, cutoffRadius, rhoInvDx, rhoSplineValues);
#endif
}

// accelerate for calculating df in computing eam potential.
inline void accelerateEamDfCalc(int *df_n, AtomElement ***atoms, double *cutoffRadius,
                                double *dfSplineInvDx, double *dfSplineValues) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, eam_df_calc)(df_n, atoms, cutoffRadius, dfSplineInvDx, dfSplineValues);
#endif
}

/**
 * accelerate for calculating force in computing eam potential.
 * // fixme many atom types.
 */
inline void accelerateEamForceCalc(int *phi_n, AtomElement ***atoms,
                                   double *cutoffRadius, double *phiSplineInvDx,
                                   double *phiSplineValues, double *rhoSplineValues) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, eam_force_calc)(phi_n, atoms, cutoffRadius, phiSplineInvDx, phiSplineValues, rhoSplineValues);
#endif
}

#endif //CRYSTAL_MD_HARDWARE_ACCELERATE_H
