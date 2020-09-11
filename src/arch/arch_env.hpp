//
// Created by genshen on 2018-3-4.
//

#ifndef MISA_MD_ARCH_H
#define MISA_MD_ARCH_H

#include <iostream>
#include "arch_building_config.h"

#include "arch_imp.h"


// initial of different platforms or different hardware architectures.
void archEnvInit() {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, env_init)();
#endif
}

void archEnvFinalize() {
#ifdef ACCELERATE_ENABLED
    ARCH_PREFIX(ARCH_NAME, env_clean)();
#endif
}

#endif //MISA_MD_ARCH_H
