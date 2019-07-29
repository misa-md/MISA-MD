//
// Created by genshen on 2019-07-29.
//

#ifndef CRYSTAL_MD_OUTPUT_INTERFACE_H
#define CRYSTAL_MD_OUTPUT_INTERFACE_H

#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"
#include "config_values.h"
#include "atom_dump.h"

/**
 * @brief abstract class of simulation step callback.
 */
class OutputInterface {
public:
    explicit OutputInterface(const Output output);

    virtual void prepareOutput(const comm::BccDomain p_domain);

    /**
     * this will be call in each output step.
     * @param time_step current time step
     * @param atom_list list of lattice atoms.
     * @param inter_atom_list list of inter atoms.
     */
    virtual void onOutputStep(const unsigned long time_step, AtomList *atom_list, InterAtomList *inter_atom_list);

    /**
     * this will be called before collision step.
     * @param time_step current collision time step
     * @param atom_list list of lattice atoms.
     * @param inter_atom_list list of inter atoms.
     */
    virtual void beforeCollision(const unsigned long time_step, AtomList *atom_list,
                                 InterAtomList *inter_atom_list);

    /**
     * this will be call when all time steps finished.
     * @param time_step current time step
     */
    virtual void onAllOut(const unsigned long time_step);

protected:
    /**
     * output configures
     */
    const Output output_config;

    /**
     * dump instance if it is not outputting by frame (only one file shared by all processors).
     */
    AtomDump *dumpInstance;

    /**
     * time of outputting
     */
    double totalDumpTime = 0;

    // atom boundary in array.
    _type_lattice_coord begin[DIMENSION];
    _type_lattice_coord end[DIMENSION];
    _type_lattice_size atoms_size;
};

#endif //CRYSTAL_MD_OUTPUT_INTERFACE_H
