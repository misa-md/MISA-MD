//
// Created by genshen on 2018-12-27.
//

#include <pack/particledata.h>
#include <pack/lat_particle_data.h>
#include "mpi_data_types.h"

MPI_Datatype mpi_types::_mpi_Particle_data;
MPI_Datatype mpi_types::_mpi_latParticle_data;

MPI_Datatype mpi_types::mpi_type_atoms_count = MPI_UNSIGNED_LONG;

void mpi_types::setInterMPIType() {
    particledata::setMPIType(_mpi_Particle_data);  // todo move code to other place?
    LatParticleData::setMPIType(_mpi_latParticle_data);
}

void mpi_types::unsetInterMPIType() {
    MPI_Type_free(&_mpi_Particle_data);
    MPI_Type_free(&_mpi_latParticle_data);
}
