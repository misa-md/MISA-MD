//
// Created by genshen on 2019-06-10.
//

#include <system_configuration.h>
#include <logs/logs.h>
#include "md_simulation.h"

MDSimulation::MDSimulation(ConfigValues *p_config_values) : simulation(p_config_values) {}

void MDSimulation::beforeStep(const unsigned long step) {
    kiwi::logs::s(MASTER_PROCESSOR, "simulation", "simulating steps: {}/{}\r",
                  _simulation_time_step + 1, pConfigVal->timeSteps);
#ifdef MD_DEV_MODE
    kiwi::logs::d("count", "real:{}--inter:{}\n", _atom->realAtoms(), _atom->getInterList()->nlocalinter);
#endif
}

void MDSimulation::postStep(const unsigned long step) {
#ifdef MD_DEV_MODE
    {
        const double e = configuration::kineticEnergy(_atom->getAtomList(), _atom->getInterList(),
                                                      configuration::ReturnMod::All, 0);
        const _type_atom_count n = 2 * pConfigVal->phaseSpace[0] *
                                   pConfigVal->phaseSpace[1] * pConfigVal->phaseSpace[2];
        const double T = configuration::temperature(e, n);
        kiwi::logs::d(MASTER_PROCESSOR, "energy", "e = {}, T = {}.\n", e, T);
    }
#endif
}
