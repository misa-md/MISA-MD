//
// Created by genshen on 5/31/18.
//

#include "atom_info_dump.h"


void atom_dump::registerAtomDumpMPIDataType() {
    /* create a type for struct car */
    const int n_items = 4;
    int blocklengths[] = {1, 1, 1, 2 * DIMENSION};
    MPI_Datatype types[] = {
            MPI_UNSIGNED_LONG, // id
            MPI_UNSIGNED_LONG, // step
            MPI_INT,  // type
            MPI_DOUBLE // location and velocity
    };
    AtomInfoDump info;
//    MPI_Aint offsets[] = {
//            offsetof(AtomInfoDump, id),
//            offsetof(AtomInfoDump, step),
//            offsetof(AtomInfoDump, type),
//            offsetof(AtomInfoDump, atom_location)
//    };
    MPI_Aint offsets[4];
    MPI_Get_address(&info.id, offsets);
    MPI_Get_address(&info.step, offsets + 1);
    MPI_Get_address(&info.type, offsets + 2);
    MPI_Get_address(&info.atom_location[0], offsets + 3);

    MPI_Type_create_struct(n_items, blocklengths, offsets, types, &atom_dump::mpi_dump_type);
    MPI_Type_commit(&atom_dump::mpi_dump_type);
}
