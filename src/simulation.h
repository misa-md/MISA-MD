//
// Created by baihe in 2016-12-22.
//

#ifndef CRYSTAL_MD_SIMULATION_H
#define CRYSTAL_MD_SIMULATION_H

#include <mpi.h>
#include <cstring>

#include "domain.h"
#include "integrator.h"
#include "create_atom.h"
#include "input.h"
#include "eam.h"
#include "toml_config.h"

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

    void prepareForStart(int rank);

    void simulate();

    void finalize();

    void initEamPotential(string file_type);

    void eamBCastPotential(int rank);

    void eamPotentialInterpolate();

    void grab(FILE *fptr, int n, double *list);

    void output();

    void exit(int exitcode);

private:
    /**
     * the time steps the program have simulated.
     */
    unsigned long _simulation_time_step;

    /**
     * pointer to config data.
     */
    ConfigValues *pConfigVal; // todo config value.
    kiwi::IOWriter *writer = nullptr; // io writer for writing a shared file using mpi-IO lib.

    Domain *_domain_decomposition; //仅rank==0的进程有效
   // GlobalDomain *_domain;  //仅rank==0的进程有效 // todo ??
    atom *_atom;
    integrator *_integrator;
    create_atom *_createatom;
    input *_input;  // 从文件读取原子坐标,速度信息
    eam *_pot;

    bool _finalCheckpoint;
};

#endif //CRYSTAL_MD_SIMULATION_H
