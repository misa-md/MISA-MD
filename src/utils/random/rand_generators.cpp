//
// Created by genshen on 2019-02-10.
//

#include "rand_generators.h"


#ifdef RAND_LCG
std::minstd_rand md_rand::rng;

void md_rand::seed(const uint32_t seed){
    rng.seed(seed);
}

uint32_t md_rand::rand32(){
     return rng();
}

#endif

#ifdef RAND_MT
std::mt19937 md_rand::rng;

void md_rand::seed(const uint32_t seed){
    rng.seed(seed);
}

uint32_t md_rand::rand32(){
     return rng();
}

#endif

#ifdef RAND_STC
std::ranlux24 md_rand::rng;

void md_rand::seed(const uint32_t seed){
    rng.seed(seed);
}

uint32_t md_rand::rand32(){
     return rng();
}

#endif

#ifdef RAND_XOSHIRO
util::random::xoroshiro128_plus md_rand::rng;

void md_rand::seed(const uint32_t seed){
     rng.init_seed(seed);
}

uint32_t md_rand::rand32(){
     return rng();
}

#endif

#ifdef RAND_LEGACY
md_rand::LegacyRand rng;

/**
 * set seed for legacy random.
 */
void md_rand::seed(const uint32_t seed) {
    rng.seed(seed);
}

uint32_t md_rand::rand32() {
    return rng.rand();
}

#endif

#ifdef RAND_LINUX_REAL
#endif