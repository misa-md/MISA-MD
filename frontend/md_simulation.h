//
// Created by genshen on 2019-06-10.
//

#ifndef CRYSTALMD_MD_SIMULATION_H
#define CRYSTALMD_MD_SIMULATION_H


#include <simulation.h>

class MDSimulation : public simulation {
public:
    explicit MDSimulation(ConfigValues *p_config_values);

    /**
     * this callback function will be called before a simulation step.
     * @param step current simulation step, starting from 0.
     */
    void beforeStep(const unsigned long step) override;

    /**
     * this callback function will be called after a simulation step.
     * @param step current simulation step, starting from 0.
     */
    void postStep(const unsigned long step) override;;

};


#endif //CRYSTALMD_MD_SIMULATION_H
