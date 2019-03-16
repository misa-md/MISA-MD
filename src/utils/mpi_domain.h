//
// Created by genshen on 9/6/18.
//

#ifndef CRYSTAL_MD_SIMULATION_DOMAIN_H
#define CRYSTAL_MD_SIMULATION_DOMAIN_H

#include <utils/mpi_utils.h>
#include <types_define.h>

/**
 * the simulation domain.
 */
struct MPIDomain {
    static kiwi::mpi_process sim_processor;

    /**
     * convert to comm::mpi_process type defined in lib comm.
     * @return copy of converted comm::mpi_process.
     */
    static comm::mpi_process toCommProcess() {
        return comm::mpi_process{
                MPIDomain::sim_processor.own_rank,
                MPIDomain::sim_processor.all_ranks,
                MPIDomain::sim_processor.comm,
        };
    }
};


#endif //CRYSTAL_MD_SIMULATION_DOMAIN_H
