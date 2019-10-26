//
// Created by genshen on 2019-07-29.
//

#ifndef CRYSTAL_MD_OUTPUT_BASE_INTERFACE_H
#define CRYSTAL_MD_OUTPUT_BASE_INTERFACE_H

#include <comm/domain/bcc_domain.h>
#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"
#include "../config_values.h"

/**
 * @brief abstract class of output atoms in simulation step.
 */
class OutputBaseInterface {
public:
    /**
     * initialize output interface by output config and domain.
     * @param config output config
     * @param domain bcc domain to specific start and end boundary of box.
     */
    explicit OutputBaseInterface(const Output config, const comm::BccDomain domain)
            : output_config(config) {
        // atom boundary in array.
        begin[0] = domain.dbx_sub_box_lattice_region.x_low - domain.dbx_ghost_ext_lattice_region.x_low;
        begin[1] = domain.dbx_sub_box_lattice_region.y_low - domain.dbx_ghost_ext_lattice_region.y_low;
        begin[2] = domain.dbx_sub_box_lattice_region.z_low - domain.dbx_ghost_ext_lattice_region.z_low;
        end[0] = begin[0] + domain.dbx_sub_box_lattice_size[0];
        end[1] = begin[1] + domain.dbx_sub_box_lattice_size[1];
        end[2] = begin[2] + domain.dbx_sub_box_lattice_size[2];
        atoms_size = domain.dbx_sub_box_lattice_size[0] *
                     domain.dbx_sub_box_lattice_size[1] *
                     domain.dbx_sub_box_lattice_size[2];
    }

    virtual ~OutputBaseInterface() {};

    virtual void prepareOutput(const comm::BccDomain p_domain) = 0;

    /**
     * this will be call in each output step.
     * @param time_step current time step
     * @param atom_list list of lattice atoms.
     * @param inter_atom_list list of inter atoms.
     */
    virtual void onOutputStep(const unsigned long time_step, AtomList *atom_list, InterAtomList *inter_atom_list) = 0;

    /**
     * this will be called before collision step.
     * @param time_step current collision time step
     * @param atom_list list of lattice atoms.
     * @param inter_atom_list list of inter atoms.
     */
    virtual void beforeCollision(const unsigned long time_step, AtomList *atom_list,
                                 InterAtomList *inter_atom_list) = 0;

    /**
     * this will be call when all time steps finished.
     * @param time_step current time step
     */
    virtual void onAllOut(const unsigned long time_step) = 0;

protected:
    /**
     * output configures
     */
    const Output output_config;

    // atom boundary in array.
    _type_lattice_coord begin[DIMENSION];
    _type_lattice_coord end[DIMENSION];
    _type_lattice_size atoms_size;
};

#endif //CRYSTAL_MD_OUTPUT_BASE_INTERFACE_H
