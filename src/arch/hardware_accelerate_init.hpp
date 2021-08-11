//
// Created by chu genshen on 2021/8/11.
//

#ifndef MISA_MD_HARDWARE_ACCELERATE_INIT_HPP
#define MISA_MD_HARDWARE_ACCELERATE_INIT_HPP

#include <comm/domain/bcc_domain.h>

#include "atom/atom_element.h"
#include "atom/neighbour_index.h"
#include "arch_building_config.h"
#include "arch_imp.h"

// callback function for hardware acceleration when domain is created.
// about const &, see: https://stackoverflow.com/questions/9637856/why-is-const-int-faster-than-const-int/9637951#9637951
inline void archAccDomainInit(const comm::BccDomain *domain) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, domain_init)(domain);
#endif
}

// callback function for acceleration, when neighbor offset indexes are created.
inline void archAccNeiOffsetInit(const NeighbourIndex<AtomElement> *nei_offset) {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, nei_offset_init)(nei_offset);
#endif
}

#endif //MISA_MD_HARDWARE_ACCELERATE_INIT_HPP
