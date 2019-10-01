//
// Created by genshen on 2019-02-10.
//

#include <random>
#include "random.h"

#ifndef MD_DEV_MODE

#include "rand_generators.h"

#endif

void md_rand::initSeed(const uint32_t seed) {
#ifndef MD_DEV_MODE
    if (seed == seed_auto) {
        std::random_device rd;
        md_rand::seed(rd());
    } else {
        md_rand::seed(seed);
    }
#else
    md_rand::seed(seed);
#endif
}

uint32_t md_rand::rand32(const uint32_t low, const uint32_t high) {
#ifdef MD_DEV_MODE
    return rand() % (high - low) + low;
#else
    return rand32() % (high - low) + low;
#endif
}

double md_rand::random() {
#ifdef MD_DEV_MODE
    return (double) rand() / RAND_MAX;
#else
    //    return (r >> 11) * (1.0 / (UINT64_C(1) << 53));
    return rand32() * (1.0 / rng.max());
#endif
}
