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

// initialize domain for hardware accelerate.
// about const &, see: https://stackoverflow.com/questions/9637856/why-is-const-int-faster-than-const-int/9637951#9637951
inline void accelerateInitDomain(const comm::BccDomain *domain) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, accelerate_init)(domain);
#endif
}

// it runs after atom and boxes creation, but before simulation running.
// which can initialize potential on device side.
inline void acceleratePotInit(eam *_pot) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, pot_init)(_pot);
#endif
}

// accelerate for calculating electron_density in computing eam potential.
inline void accelerateEamRhoCalc(eam *pot, AtomElement *atoms, const double cutoff_radius) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, eam_rho_calc)(pot, atoms, cutoff_radius);
#endif
}

// accelerate for calculating df in computing eam potential.
inline void accelerateEamDfCalc(eam *pot, AtomElement *atoms, const double cutoff_radius) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, eam_df_calc)(pot, atoms, cutoff_radius);
#endif
}

/**
 * accelerate for calculating force in computing eam potential.
 */
inline void accelerateEamForceCalc(eam *pot, AtomElement *atoms, const double cutoff_radius) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, eam_force_calc)(pot, atoms, cutoff_radius);
#endif
}

#endif //CRYSTAL_MD_HARDWARE_ACCELERATE_H
