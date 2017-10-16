#include <mpi.h>
#include "latparticledata.h"

void latparticledata::setMPIType(MPI_Datatype &sendPartType) {
	int blocklengths[] = { 1, 3 }; // 1 int value (type), 3 double values (3r)
	MPI_Datatype types[] = { MPI_INT, MPI_DOUBLE };

	MPI_Aint displacements[2];
	latparticledata pdata_dummy;
#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
	MPI_Get_address(&pdata_dummy, displacements);
	MPI_Get_address(&pdata_dummy.r[0], displacements + 1);
#else
	MPI_Address(&pdata_dummy, displacements);
	MPI_Address(&pdata_dummy.r[0], displacements + 1);
#endif
	MPI_Aint base = displacements[0];
	for (int i = 0; i < 2; i++)
		displacements[i] -= base;

#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
	MPI_Type_create_struct(2, blocklengths, displacements, types, &sendPartType);
#else
    MPI_Type_struct(2, blocklengths, displacements, types, &sendPartType);
#endif
	MPI_Type_commit(&sendPartType);
}

latparticledata::latparticledata() : type(-1) {
	for (int i = 0; i < 3; i++ ) {
		r[i] = 0.0;
	}
}