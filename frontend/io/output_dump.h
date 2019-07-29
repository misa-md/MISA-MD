//
// Created by genshen on 2019-07-29.
//

#ifndef CRYSTALMD_OUTPUT_DUMP_H
#define CRYSTALMD_OUTPUT_DUMP_H

#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"
#include "config_values.h"
#include "output_base_interface.h"

/**
 * @brief abstract class of simulation step callback.
 */
class OutputDump : public OutputBaseInterface {
public:
    explicit OutputDump(const Output output, const comm::BccDomain p_domain);

    void prepareOutput(const comm::BccDomain p_domain) override;

    /**
     * this will be call in each output step.
     * @param time_step current time step
     * @param atom_list list of lattice atoms.
     * @param inter_atom_list list of inter atoms.
     */
    void onOutputStep(const unsigned long time_step, AtomList *atom_list,
                      InterAtomList *inter_atom_list) override;

    /**
     * this will be called before collision step.
     * @param time_step current collision time step
     * @param atom_list list of lattice atoms.
     * @param inter_atom_list list of inter atoms.
     */
    void beforeCollision(const unsigned long time_step, AtomList *atom_list,
                         InterAtomList *inter_atom_list) override;

    /**
     * this will be call when all time steps finished.
     * @param time_step current time step
     */
    void onAllOut(const unsigned long time_step) override;

protected:
    /**
     * time of outputting
     */
    double total_dump_time = 0;

    static void dump(const unsigned long time_step, AtomList *atom_list, InterAtomList *inter_atom_list);
};

#endif //CRYSTALMD_OUTPUT_DUMP_H
