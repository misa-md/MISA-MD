#include "particledata.h"

void particledata::setMPIType(MPI_Datatype &sendPartType) {
    int blocklengths[] = {1, 1, 6}; // 1 unsLong value (id), 1 int value (cid), 6 double values (3r, 3v)
    MPI_Datatype types[] = {MPI_UNSIGNED_LONG, MPI_INT, MPI_DOUBLE};

    MPI_Aint displacements[3];
    particledata pdata_dummy;
#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
    MPI_Get_address(&pdata_dummy, displacements);
    MPI_Get_address(&pdata_dummy.type, displacements + 1);
    MPI_Get_address(&pdata_dummy.r[0], displacements + 2);
#else
    MPI_Address(&pdata_dummy, displacements);
    MPI_Address(&pdata_dummy.type, displacements + 1);
    MPI_Address(&pdata_dummy.r[0], displacements + 2);
#endif
    MPI_Aint base = displacements[0];
    for (int i = 0; i < 3; i++) {
        displacements[i] -= base;
    }

#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
    MPI_Type_create_struct(3, blocklengths, displacements, types, &sendPartType);
#else
    MPI_Type_struct(3, blocklengths, displacements, types, &sendPartType);
#endif
    MPI_Type_commit(&sendPartType);
}

particledata::particledata() : id(0), type(atom_type::Fe) {
    for (int i = 0; i < 3; i++) {
        r[i] = 0.0;
        v[i] = 0.0;
    }
}
