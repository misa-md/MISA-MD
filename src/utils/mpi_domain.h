//
// Created by genshen on 9/6/18.
//

#ifndef CRYSTAL_MD_SIMULATION_DOMAIN_H
#define CRYSTAL_MD_SIMULATION_DOMAIN_H

#include <utils/mpi_utils.h>

/**
 * the simulation domain.
 */
struct MPIDomain {
    static kiwi::mpi_process sim_processor;
};


#endif //CRYSTAL_MD_SIMULATION_DOMAIN_H
