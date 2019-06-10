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

#ifdef MD_DEV_MODE
    {
        unsigned long count[2] = {0, 0};
        count[0] = _atom->realAtoms();
        count[1] = _atom->getInterList()->nlocalinter;
        unsigned long global_count[2] = {0, 0};
        unsigned long &real_atoms = global_count[0], &inter_atoms = global_count[1];
        MPI_Reduce(count, global_count, 2, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_PROCESSOR, MPI_COMM_WORLD);

        kiwi::logs::d("count", "real:{}--inter:{}\n", count[0], count[1]);
        kiwi::logs::d(MASTER_PROCESSOR, "count", "global_real:{}--global_inter:{}\n", real_atoms, inter_atoms);
    }
#endif
}

void MDSimulation::onForceSolved(const unsigned long step) {
#ifdef MD_DEV_MODE
    forceChecking();
#endif
}

void MDSimulation::print_force(const std::string filename, int step) {
    std::ofstream outfile;
    outfile.open(filename);

    _atom->atom_list->foreachSubBoxAtom(
            [&outfile](AtomElement &_atom_ref) {
                outfile << _atom_ref.id << "\t"
                        << _atom_ref.x[0] << "\t" << _atom_ref.x[1] << "\t" << _atom_ref.x[2] << "\t"
                        << _atom_ref.v[0] << "\t" << _atom_ref.v[1] << "\t" << _atom_ref.v[2] << "\t"
                        << _atom_ref.f[0] << "\t" << _atom_ref.f[1] << "\t" << _atom_ref.f[2] << std::endl;
            }
    );
    outfile << "inter atoms" << std::endl;
    for (_type_inter_list::iterator inter_it = _atom->inter_atom_list->inter_list.begin();
         inter_it != _atom->inter_atom_list->inter_list.end(); inter_it++) {
        AtomElement &_atom_ref = *inter_it;
        outfile << std::endl << _atom_ref.id << "\t"
                << _atom_ref.x[0] << "\t" << _atom_ref.x[1] << "\t" << _atom_ref.x[2] << "\t"
                << _atom_ref.v[0] << "\t" << _atom_ref.v[1] << "\t" << _atom_ref.v[2] << "\t"
                << _atom_ref.f[0] << "\t" << _atom_ref.f[1] << "\t" << _atom_ref.f[2] << std::endl;
    }
    outfile.close();
}


#ifdef MD_DEV_MODE

void MDSimulation::forceChecking() {
    auto forces = configuration::systemForce(_atom->getAtomList(), _atom->getInterList());
    double fx[3] = {forces[0], forces[1], forces[2]};
    double fx_2[3] = {0.0, 0.0, 0.0};
    MPI_Allreduce(fx, fx_2, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    kiwi::logs::d("force2", "({}, {}, {})\n\n", forces[0], forces[1], forces[2]);
    kiwi::logs::d(MASTER_PROCESSOR, "sum force:", "({}, {}, {})\n\n", fx_2[0], fx_2[1], fx_2[2]);
    if (std::abs(fx_2[0]) > 0.00001) {
        char filename[20];
        sprintf(filename, "force_%d.txt", kiwi::mpiUtils::global_process.own_rank);
        print_force(filename, _simulation_time_step);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
}

#endif
