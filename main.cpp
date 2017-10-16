#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include "simulation.h"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int world_size = 1;
    int ownrank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &ownrank);
    if (ownrank == 0)
        std::cout << "Running with " << world_size << " MPI processes." << std::endl;
    simulation simulation;

    //读取数据或者创建数据
    simulation.createboxandatom();

    //准备模拟
    simulation.prepare_start(ownrank);

    if (ownrank == 0)
        std::cout << "Initializing simulation" << std::endl;
    //开始模拟
    simulation.simulate();
    if (ownrank == 0)
        std::cout << "finalizing simulation" << std::endl;
    //模拟结束
    simulation.finalize();

    MPI_Finalize();
    return 0;
}
