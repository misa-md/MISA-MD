//
// Created by baihe back to 2015-12-26.
//

#ifndef PARTICLEDATA_H_
#define PARTICLEDATA_H_

#include <mpi.h>
#include "../atom/atom_types.h"

class particledata {
public:
    // 定义一个数据类型，用于MPI数据传输
    static void setMPIType(MPI_Datatype &sendPartType);

    particledata();

    unsigned long id;
    atom_type::atom_type type;
    double r[3];
    double v[3];
};

#endif /* PARTICLEDATA_H_ */
