//
// Created by genshen on 2018-3-4.
//

#ifndef CRYSTALMD_ARCH_H
#define CRYSTALMD_ARCH_H

#include <iostream>
#include "pre_config.h"

#ifdef ARCH_SUNWAY

#include "arch/sunway/sunway.h"

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

#endif //CRYSTALMD_ARCH_H
