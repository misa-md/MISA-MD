//
// Created by genshen on 2018-12-27.
//

#ifndef MISA_MD_MPI_DATA_TYPES_H
#define MISA_MD_MPI_DATA_TYPES_H


#include <mpi.h>

namespace mpi_types {

    /** mpi data structure for communicating inter atoms. **/
    extern MPI_Datatype _mpi_Particle_data;
    extern MPI_Datatype _mpi_latParticle_data;

    // todo call me in initialization
    void setInterMPIType();

    // todo call me when finished
    void unsetInterMPIType();
};


#endif //MISA_MD_MPI_DATA_TYPES_H
