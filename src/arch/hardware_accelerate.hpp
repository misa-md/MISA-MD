//
// Created by genshen on 2018-3-5.
//

#ifndef MISA_MD_HARDWARE_ACCELERATE_H
#define MISA_MD_HARDWARE_ACCELERATE_H

#include <eam.h>
#include <comm/domain/bcc_domain.h>
#include <args.hpp>

#include "atom/atom_element.h"
#include "arch_building_config.h"
#include "arch_imp.h"

// check whether it has accelerate hardware to be used, for example GPU, MIC(Xeon Phi), or sunway slave cores.
inline bool isArchAccSupport() {
#ifdef ACCELERATE_ENABLED
    return true; // sunway and other hardware.
#else
    return false;
#endif
}

// set command line interface options.
inline void archCliOptions(args::ArgumentParser &parser) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, cli_options)(parser);
#endif
}

inline bool archCliOptionsParse(args::ArgumentParser &parser) {
#ifdef ACCELERATE_ENABLED
    return ARCH_PREFIX(ARCH_NAME, cli_options_parse)(parser);
#else
    return true;
#endif
}

// api for creating lattice atoms memory.
// It will recreate using `malloc/new` after this call if it returns nullptr.
template<typename T>
inline T *archCreateAtomsMemory(_type_atom_count size_x, _type_atom_count size_y, _type_atom_count size_z) {
#ifdef ACCELERATE_ENABLED
    return ARCH_PREFIX(ARCH_NAME, create_atoms_mem)(size_x, size_y, size_z);
#else
    return nullptr;
#endif
}

// release created lattice atoms.
template<typename T>
inline bool archReleaseAtomsMemory(T *atoms) {
#ifdef ACCELERATE_ENABLED
    return ARCH_PREFIX(ARCH_NAME, release_atoms_mem)(atoms);
#else
    return false;
#endif
}

// it runs after atom and boxes creation, but before simulation running.
// which can initialize potential on device side.
inline void archAccPotInit(eam *_pot) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, pot_init)(_pot);
#endif
}

// accelerate for calculating electron_density in computing eam potential.
inline void archAccEamRhoCalc(eam *pot, AtomElement *atoms, const double cutoff_radius) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, eam_rho_calc)(pot, atoms, cutoff_radius);
#endif
}

// accelerate for calculating df in computing eam potential.
inline void archAccEamDfCalc(eam *pot, AtomElement *atoms, const double cutoff_radius) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, eam_df_calc)(pot, atoms, cutoff_radius);
#endif
}

/**
 * accelerate for calculating force in computing eam potential.
 */
inline void archAccEamForceCalc(eam *pot, AtomElement *atoms, const double cutoff_radius) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, eam_force_calc)(pot, atoms, cutoff_radius);
#endif
}

#endif //MISA_MD_HARDWARE_ACCELERATE_H
