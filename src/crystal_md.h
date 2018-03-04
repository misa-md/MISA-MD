//
// Created by gensh(genshenchu@gmail.com) on 2017/10/19.
//

#ifndef CRYSTALMD_CRYSTAL_MD_H
#define CRYSTALMD_CRYSTAL_MD_H

#include "simulation.h"
#include "config.h"

class crystalMD {

public:

    crystalMD(int argc, char **argv);

    // initial MPI env; read terminal args and config file
    //初始化MPI环境,读取命令行参数,然后根据参数读取配置文件并解析config.
    bool initialize();

    //create boxes and atoms for later simulation.
    //进行区域分解,创建原子
    bool prepare();

    //运行模拟
    void run();

    void destroy();

    void detach();

private:
    int argc = 0;
    char **argv;

    short argvStatus = 0;
    config *cp;
    simulation *simu;

    bool runtimeEnvInitialize();
};

#endif //CRYSTALMD_CRYSTAL_MD_H
