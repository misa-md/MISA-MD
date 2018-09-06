//
// Created by genshen on 2018-05-19.
//

#include <gtest/gtest.h>
#include <potential/eam.h>
#include <utils/mpi_utils.h>
#include "utils/mpi_domain.h"

TEST(phi_test, phi_sync) {
    eam *_eam = new eam();
    if (MPIDomain::sim_processor.own_rank == MASTER_PROCESSOR) {
        _eam->initElementN(1);
        double data[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
        _eam->eam_phi.append(atom_type::Fe, atom_type::Fe, 6, 0, 1.0, data);
        _eam->electron_density.append(atom_type::Fe, 6, 0, 1.0, data);
        // only the data is added, initInterpolationObject can be called. todo so some checks must be added.
        _eam->embedded.append(atom_type::Fe, 6, 0, 1.0, data);
    }
    _eam->eamBCast(MPIDomain::sim_processor.own_rank); // sync to other processors.
    _eam->interpolateFile();
    delete _eam;
}
