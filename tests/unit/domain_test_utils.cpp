//
// Created by genshen on 5/12/18.
//

#include "domain_test_utils.h"

Domain *getDomainInstance(int64_t space[3], double lattice_const, double cutoff_radius) {
    return Domain::Builder()
            .setPhaseSpace(space)
            .setCutoffRadius(cutoff_radius)
            .setLatticeConst(lattice_const)
            .build();
}
