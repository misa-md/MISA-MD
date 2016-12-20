#ifndef LATPARTICLEDATA_H_
#define LATPARTICLEDATA_H_

#include <mpi.h>

class latparticledata {
public:
	// 定义一个数据类型，用于MPI数据传输
	static void setMPIType(MPI_Datatype &sendPartType);

	latparticledata();

	int type;
	double r[3];
};

#endif /* LATPARTICLEDATA_H_ */