//
// Created by genshen on 2019-02-10.
//

#include "rand_generators.h"

md_rand::type_rng md_rand::rng;

// all rand number generation engines share the same implementation to generate a new random number.
uint32_t md_rand::rand32() {
    return rng();
}

void md_rand::seed(const uint32_t seed) {
    rng.seed(seed);
}

#ifdef RAND_LINUX_REAL
#endif
