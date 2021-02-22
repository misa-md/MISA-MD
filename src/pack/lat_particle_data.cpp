#include "lat_particle_data.h"

void LatParticleData::setMPIType(MPI_Datatype &sendPartType) {
#ifdef DEV_MD_COMM_INC_ATOM_ID
    constexpr int data_types = 3;
    int blocklengths[] = {1, 3, 1}; // 1 int value (type), 3 double values (3r), 1 long type (id)
    MPI_Datatype types[] = {MPI_INT, MPI_DOUBLE, MPI_UNSIGNED_LONG};
#else
    constexpr int data_types = 2;
    int blocklengths[] = {1, 3}; // 1 int value (type), 3 double values (3r)
    MPI_Datatype types[] = {MPI_INT, MPI_DOUBLE};
#endif

    MPI_Aint displacements[data_types];
    LatParticleData pdata_dummy;

#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
    MPI_Get_address(&pdata_dummy, displacements);
    MPI_Get_address(&pdata_dummy.r[0], displacements + 1);
#ifdef DEV_MD_COMM_INC_ATOM_ID
    MPI_Get_address(&pdata_dummy.id, displacements + 2);
#endif
#else
    MPI_Address(&pdata_dummy, displacements);
    MPI_Address(&pdata_dummy.r[0], displacements + 1);
#ifdef DEV_MD_COMM_INC_ATOM_ID
    MPI_Address(&pdata_dummy.id, displacements + 2);
#endif
#endif

    MPI_Aint base = displacements[0];
    for (int i = 0; i < data_types; i++) {
        displacements[i] -= base;
    }

#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
    MPI_Type_create_struct(data_types, blocklengths, displacements, types, &sendPartType);
#else
    MPI_Type_struct(data_types, blocklengths, displacements, types, &sendPartType);
#endif
    MPI_Type_commit(&sendPartType);
}

LatParticleData::LatParticleData() : type(atom_type::Fe) {
    for (int i = 0; i < 3; i++) {
        r[i] = 0.0;
    }
}
