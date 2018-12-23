//
// Created by baihe in 2016-12-22.
//

#ifndef CRYSTAL_MD_SIMULATION_H
#define CRYSTAL_MD_SIMULATION_H

#include <mpi.h>
#include <cstring>
#include <io/io_writer.h>
#include <eam.h>

#include "toml_config.h"
#include "newton_motion.h"
#include "input.h"
#include "atom_dump.h"

class simulation {
public:

    simulation();

    ~simulation();

    /**
     * Denote N as the count of all processors.
     * {@memberof domainDecomposition} will divide the simulation box into N parts,
     * we call each part as a sub-box.
     * And each sub-box will bind to a processor.
     */
    void createDomainDecomposition();

    void createAtoms();

    void prepareForStart();

    void simulate();

    void finalize();

    /**
     * start to dump atoms to file
     * @param time_step current time step
     * @param before_collision true for dumping atoms before collision
     */
    void output(size_t time_step, bool before_collision = false);

    void abort(int exitcode);

private:
    /**
     * the time steps the program have simulated.
     */
    unsigned long _simulation_time_step;

    /**
     * pointer to config data.
     */
    ConfigValues *pConfigVal;
    Domain *_p_domain; //仅rank==0的进程有效
    // GlobalDomain *p_domain;  //仅rank==0的进程有效 // todo ??
    atom *_atom;
    NewtonMotion *_newton_motion;

    input *_input;  // 从文件读取原子坐标,速度信息
    eam *_pot; // eam potential

    bool _finalCheckpoint;
};

#endif //CRYSTAL_MD_SIMULATION_H
