//
// Created by genshen on 5/12/18.
//

#include "domain_test_utils.h"

Domain *getDomainInstance(int64_t space[3], double lattice_const, double cutoff_radius) {
    return (new Domain(space, lattice_const, cutoff_radius))
            ->decomposition()
            ->createGlobalDomain()
            ->createSubBoxDomain();
}
