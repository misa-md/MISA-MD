//
// Created by genshen on 2020/5/22.
//

#ifndef CRYSTAL_MD_ARCH_MACROS_H
#define CRYSTAL_MD_ARCH_MACROS_H

#include "arch_building_config.h"

#define _ARCH_PREFIX_(arch_name, func_name)   arch_name ## _ ## func_name
#define ARCH_PREFIX(arch_name, func_name) _ARCH_PREFIX_(arch_name, func_name)

#endif //CRYSTAL_MD_ARCH_MACROS_H
