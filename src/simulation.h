//
// Created by baihe in 2016-12-22.
//

#ifndef CRYSTAL_MD_SIMULATION_H
#define CRYSTAL_MD_SIMULATION_H

#include <mpi.h>
#include <cstring>
#include <io/io_writer.h>
#include <eam.h>

#include "newton_motion.h"
#include "input.h"

class simulation {
public:

    simulation();

    virtual ~simulation();

    /**
     * Denote N as the count of all processors.
     * {@memberof domainDecomposition} will divide the simulation box into N parts,
     * we call each part as a sub-box.
     * And each sub-box will bind to a processor.
     * @param phase_space the lattice count in each dimension.
     * @param lattice_const lattice constance
     * @param cutoff_radius_factor cutoff radius factor = cutoff/lattice_const
     */
    void createDomain(const int64_t phase_space[DIMENSION],
                      const double lattice_const, const double cutoff_radius_factor);

    /**
     * initialize atoms position and velocity by given arguments.
     * @param phase_space the lattice count in each dimension.
     * @param lattice_const lattice const.
     * @param init_step_len time step length as initialized value.
     * @param create_mode whether create atom randomly.
     * @param create_seed seed to create atoms.
     * @param t_set initial temperature.
     * @param alloy_ratio ratio of alloy for each types of material.
     */
    void createAtoms(const int64_t phase_space[DIMENSION], const double lattice_const, const double init_step_len,
                     const bool create_mode, const unsigned long create_seed, const double t_set,
                     const int alloy_ratio[atom_type::num_atom_types]);

    /**
     * initialize potential function and perform the first simulation step.
     * @param pot_file_path path of potential file.
     */
    void prepareForStart(const std::string pot_file_path);

    /**
     *
     * @param steps total simulation steps.
     * @param coll_step step to perform collision.
     * @param coll_lat lattice position of PKA.
     * @param coll_dir direction of collision.
     * @param coll_pka_energy pka energy, unit eV.
     */
    void simulate(const unsigned long steps, unsigned long coll_step,
                  const _type_lattice_coord coll_lat[DIMENSION + 1], const double coll_dir[DIMENSION],
                  const double coll_pka_energy);

    void finalize();

    /**
     * this function will be called before simulation loop.
     */
    virtual void onSimulationStarted() {};

    /**
     * this function will be called after simulation loop finished.
     * @param step total simulation step.
     */
    virtual void onSimulationDone(const unsigned long step) {};

    /**
     * this callback function will be called before a simulation step.
     * @param step current simulation step, starting from 0.
     */
    virtual void beforeStep(const unsigned long step) {};

    /**
     * this callback function will be called after a simulation step.
     * @param step current simulation step, starting from 0.
     */
    virtual void postStep(const unsigned long step) {};

    /**
     * this callback function will be called after atoms' force are computed.
     * @param step current simulation step, starting from 0.
     */
    virtual void onForceSolved(const unsigned long step) {};

    void abort(int exitcode);

protected:
    /**
     * the time steps the program have simulated.
     */
    unsigned long _simulation_time_step;

    comm::BccDomain *_p_domain;
    // GlobalDomain *p_domain;  //仅rank==0的进程有效 // todo ??
    atom *_atom;
    NewtonMotion *_newton_motion;

    input *_input;  // 从文件读取原子坐标,速度信息
    eam *_pot; // eam potential

};

#endif //CRYSTAL_MD_SIMULATION_H
