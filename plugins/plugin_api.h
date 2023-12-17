//
// Created by genshen on 2023/12/13.
//

#ifndef MISA_MD_PLUGIN_API_H
#define MISA_MD_PLUGIN_API_H

#include "atom/atom_set.h"
#include "lattice/ws_utils.h"

namespace plugins {
    class IOPlugin {
    public:
        virtual ~IOPlugin() = default;

        virtual void init() = 0;

        virtual bool filter_atom(const _type_atom_location pos[DIMENSION]) = 0;
    };
}

#endif //MISA_MD_PLUGIN_API_H
