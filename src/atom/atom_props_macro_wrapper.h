//
// Created by chu genshen on 2021/8/25.
//

#ifndef MISA_MD_ATOM_PROPS_MACRO_WRAPPER_H
#define MISA_MD_ATOM_PROPS_MACRO_WRAPPER_H

#include "types/atom_types.h"
#include "md_building_config.h"

#ifdef MD_ATOM_HASH_ARRAY_MEMORY_LAYOUT_AOS
#include "atom_props_macro_aos.h"
#endif

#ifndef MD_ATOM_HASH_ARRAY_MEMORY_LAYOUT_AOS
#include "atom_props_macro_soa.h"
#endif

#endif //MISA_MD_ATOM_PROPS_MACRO_WRAPPER_H
