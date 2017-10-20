#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <mpi.h>
#include <cstring>

#include "domaindecomposition.h"
#include "atom.h"
#include "integrator.h"
#include "createatom.h"
#include "input.h"
#include "eam.h"
#include "config.h"


class simulation {
public:

    simulation();

    ~simulation();

    void domainDecomposition();

    void createBoxedAndAtoms();

    void prepareForStart(int rank);

    void simulate();

    void finalize();

    void initEamPotential(string file_type);

    void eamBcastPotential(int rank);

    void eamPotentialInterpolate();

    void grab(FILE *fptr, int n, double *list);

    void output();

    void exit(int exitcode);

private:
//    unsigned long _numberOfTimesteps;   /**< 时间步 */
    unsigned long _simstep;             /**< 程序运行实际时间步 */
//    int ownrank;

//    unsigned long collision_step;
//    int lat[4];
//    double collision_v[3];

//    string filetype;
//    string filename;

    config *cp;

    domaindecomposition *_domaindecomposition; //仅rank==0的进程有效
    domain *_domain;  //仅rank==0的进程有效
    atom *_atom;
    integrator *_integrator;
    createatom *_createatom;
    input *_input;  // 从文件读取原子坐标,速度信息
    eam *_pot;

    bool _finalCheckpoint;
};

#endif /*SIMULATION_H_*/