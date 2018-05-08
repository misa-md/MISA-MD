//
// Created by genshen on 2018-3-5.
//

#ifndef CRYSTAL_MD_HARDWARE_ACCELERATE_H
#define CRYSTAL_MD_HARDWARE_ACCELERATE_H

#include "pre_define.h"

#ifdef ARCH_SUNWAY

#include "arch/sunway/athread_accelerate.h" // sunway athread

#endif

// check whether it has accelerate hardware to be used, for example GPU, MIC(Xeon Phi), or sunway slave cores.
inline bool isAccelerateSupport() {
#ifdef ARCH_SUNWAY
    return true; // sunway.
#else // todo other hardware.
    return false;
#endif
}

// initial for hardware accelerate.
// about const &, see: https://stackoverflow.com/questions/9637856/why-is-const-int-faster-than-const-int/9637951#9637951
inline void accelerateInit(const int lolocalx, const int lolocaly, const int lolocalz,
                           const int nlocalx, const int nlocaly, const int nlocalz,
                           const int loghostx, const int loghosty, const int loghostz,
                           const int nghostx, const int nghosty, const int nghostz) {
#ifdef ARCH_SUNWAY
    athreadAccelerateInit(lolocalx, lolocaly, lolocalz, nlocalx, nlocaly, nlocalz,
                          loghostx, loghosty, loghostz, nghostx, nghosty, nghostz);
#endif
}

// it runs after atom and boxes creation, but before simulation running.
inline void beforeAccelerateRun(eam *_pot) {
#ifdef ARCH_SUNWAY
    initSpline(_pot->rho->spline, _pot->f->spline, _pot->phi->spline);
#endif
}

// accelerate for calculating rho in computing eam potential.
inline void accelerateEamRhoCalc(int *rho_n, AtomList *atom_list, double *cutoffRadius,
                                 double *rhoInvDx, double *rhoSplineValues) {
#ifdef ARCH_SUNWAY
    athreadAccelerateEamRhoCalc(rho_n, x, rho, cutoffRadius, rhoInvDx, rhoSplineValues);
#endif
}

// accelerate for calculating df in computing eam potential.
inline void accelerateEamDfCalc(int *df_n, AtomList *atom_list, double *cutoffRadius,
                                double *dfSplineInvDx, double *dfSplineValues) {
#ifdef ARCH_SUNWAY
    athreadAccelerateEamDfCalc(df_n, rho, df, cutoffRadius, dfSplineInvDx, dfSplineValues);
#endif
}

// accelerate for calculating force in computing eam potential.
inline void accelerateEamForceCalc(int *phi_n, AtomList *atom_list,
                                   double *cutoffRadius, double *phiSplineInvDx,
                                   double *phiSplineValues, double *rhoSplineValues) {
#ifdef ARCH_SUNWAY
    athreadAccelerateEamForceCalc(phi_n, x, f, df, cutoffRadius, phiSplineInvDx, phiSplineValues, rhoSplineValues);
#endif
}

#endif //CRYSTAL_MD_HARDWARE_ACCELERATE_H
