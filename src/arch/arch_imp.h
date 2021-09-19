//
// Created by genshen on 2019/11/13.
//

#ifndef MISA_MD_ARCH_IMP_H
#define MISA_MD_ARCH_IMP_H

#include "arch_building_config.h"

#ifdef ACCELERATE_ENABLED
#include <comm/domain/bcc_domain.h>
#include <eam.h>
#include <args.hpp>

#include "arch_macros.h"
#include "arch_atom_list_collection.h"
#include "atom/atom_element.h"
#include "atom/neighbour_index.h"

void ARCH_PREFIX(ARCH_NAME, cli_options)(args::ArgumentParser &parser);

// if you return false, it will abort argument parsing.
// You may print error messages in the parsing implementation.
bool ARCH_PREFIX(ARCH_NAME, cli_options_parse)(args::ArgumentParser &parser);

void ARCH_PREFIX(ARCH_NAME, env_init)();

void ARCH_PREFIX(ARCH_NAME, env_clean)();

// create atoms memory for saving lattice atoms.
// If it returns false, it means creations is failed or rejected.
bool ARCH_PREFIX(ARCH_NAME, create_atoms_mem)(void **atoms, std::size_t element_size, _type_atom_count size_x,
                                              _type_atom_count size_y, _type_atom_count size_z);

// release atoms memory of saving lattice atoms.
// If it returns false, it means destroying is failed or rejected.
bool ARCH_PREFIX(ARCH_NAME, release_atoms_mem)(void *atoms);

void ARCH_PREFIX(ARCH_NAME, domain_init)(const comm::BccDomain *domain);

void ARCH_PREFIX(ARCH_NAME, nei_offset_init)(const NeighbourIndex<_type_neighbour_index_ele> *nei_offset);

void ARCH_PREFIX(ARCH_NAME, pot_init)(eam *_pot);

void ARCH_PREFIX(ARCH_NAME, eam_rho_calc)(eam *pot, _type_atom_list_collection atoms, const double cutoff_radius);

void ARCH_PREFIX(ARCH_NAME, eam_df_calc)(eam *pot, _type_atom_list_collection atoms, const double cutoff_radius);

void ARCH_PREFIX(ARCH_NAME, eam_force_calc)(eam *pot, _type_atom_list_collection atoms, const double cutoff_radius);

#endif //ACCELERATE_ENABLED

#endif //MISA_MD_ARCH_IMP_H
