//
// Created by genshen on 2018-3-4.
//

#ifndef CRYSTAL_MD_ARCH_H
#define CRYSTAL_MD_ARCH_H

#include <iostream>
#include "pre_define.h"

#ifdef ARCH_SUNWAY

#include "arch/sunway/sunway_env.h"

#endif

// initial of different platforms or different hardware architectures.
void archEnvInit() {
#ifdef ARCH_SUNWAY
    sunwayAThreadInit();
#endif
}

void archEnvFinalize() {
#ifdef ARCH_SUNWAY
    sunwayAThreadClean();
#endif
}

#endif //CRYSTAL_MD_ARCH_H
