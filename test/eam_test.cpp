//
// Created by genshen on 2018-05-19.
//

#include <gtest/gtest.h>
#include <potential/eam.h>
#include <utils/mpi_utils.h>

TEST(phi_test, phi_sync) {
    eam *_eam = new eam();
    if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
        _eam->initElementN(3);
        double data[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
        _eam->eam_phi.append(atom_type::Fe, atom_type::Fe, 6, 0, 1.0, data);
    }
    _eam->eamBCast(kiwi::mpiUtils::own_rank); // sync to other processors.
    _eam->interpolateFile(); // fixme segment fault.
    delete _eam;
}
