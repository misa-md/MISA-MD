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
    void postStep(const unsigned long step) override;

    /**
     * this callback function will be called after atoms' force are computed.
     * @param step current simulation step, starting from 0.
     */
    void onForceSolved(const unsigned long step) override;

private:
#ifdef MD_DEV_MODE

    /**
     * check force, only works in dev mode.
     * if the sum of all atoms not 0, program will exists with error message.
     */
    void forceChecking();

#endif

};


#endif //CRYSTALMD_MD_SIMULATION_H
