//
// Created by gensh(genshenchu@gmail.com)  on 2017/4/15.
//

#include "crystal_md.h"

int main(int argc, char **argv) {
    // app's lifecycle here.
    auto *app = new crystalMD(argc, argv);
    if ((app->initialize())) {
        if (app->prepare()) {
            app->run();
            app->destroy();
        }
    }
    app->detach();
    return 0;
}
/*
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
    simulation.createBoxedAndAtoms();

    //准备模拟
    simulation.prepareForStart(ownrank);

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
*/