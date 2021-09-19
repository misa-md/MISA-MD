//
// Created by chu genshen on 2021/8/26.
//

#ifndef MISA_MD_ARCH_ATOM_LIST_COLLECTION_H
#define MISA_MD_ARCH_ATOM_LIST_COLLECTION_H

#include "md_building_config.h"
#include "atom/atom_element.h"

#ifdef MD_ATOM_HASH_ARRAY_MEMORY_LAYOUT_AOS

typedef AtomElement _type_neighbour_index_ele;

// pointer to data in hash array of atom list.
typedef struct {
    AtomElement *atoms;
} _type_atom_list_collection;

#endif //MD_ATOM_HASH_ARRAY_MEMORY_LAYOUT_AOS

#ifndef MD_ATOM_HASH_ARRAY_MEMORY_LAYOUT_AOS

typedef _type_atom_id _type_neighbour_index_ele;

// pointer to all data in hash array of atom list.
typedef struct {
    _type_atom_id *atom_ids;
    _type_atom_type_enum *atom_types;
    _type_atom_location (*atom_x)[DIMENSION];
    _type_atom_velocity (*atom_v)[DIMENSION];
    _type_atom_force (*atom_f)[DIMENSION];
    _type_atom_rho *atom_rho;
    _type_atom_rho *atom_df;
} _type_atom_list_collection;

#endif //MD_ATOM_HASH_ARRAY_MEMORY_LAYOUT_AOS

#endif //MISA_MD_ARCH_ATOM_LIST_COLLECTION_H
