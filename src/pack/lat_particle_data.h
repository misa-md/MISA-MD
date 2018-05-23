//
// Created by baihe back to 2016-06-28.
//

#ifndef CRYSTAL_MD_LAT_PARTICLE_DATA_H
#define CRYSTAL_MD_LAT_PARTICLE_DATA_H

#include <mpi.h>
#include "../atom/atom_types.h"

class LatParticleData {
public:
    // 定义一个数据类型，用于MPI数据传输
    static void setMPIType(MPI_Datatype &sendPartType);

    LatParticleData();

    atom_type::atom_type type;
    double r[3];
};

#endif //CRYSTAL_MD_LAT_PARTICLE_DATA_H
