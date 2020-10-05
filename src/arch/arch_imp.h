//
// Created by genshen on 2019/11/13.
//

#ifndef MISA_MD_ARCH_IMP_H
#define MISA_MD_ARCH_IMP_H

#include "arch_building_config.h"

#ifdef ACCELERATE_ENABLED
#include <comm/domain/bcc_domain.h>
#include <eam.h>

#include "atom/atom_element.h"
#include "atom/neighbour_index.h"

void ARCH_PREFIX(ARCH_NAME, env_init)();

void ARCH_PREFIX(ARCH_NAME, env_clean)();

void ARCH_PREFIX(ARCH_NAME, domain_init)(const comm::BccDomain *domain);

void ARCH_PREFIX(ARCH_NAME, nei_offset_init)(const NeighbourIndex<AtomElement> *nei_offset);

void ARCH_PREFIX(ARCH_NAME, pot_init)(eam *_pot);

void ARCH_PREFIX(ARCH_NAME, eam_rho_calc)(eam *pot, AtomElement *atoms, const double cutoff_radius);

void ARCH_PREFIX(ARCH_NAME, eam_df_calc)(eam *pot, AtomElement *atoms, const double cutoff_radius);

void ARCH_PREFIX(ARCH_NAME, eam_force_calc)(eam *pot, AtomElement *atoms, const double cutoff_radius);

#endif //ACCELERATE_ENABLED

#endif //MISA_MD_ARCH_IMP_H
