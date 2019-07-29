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

    if (step == pConfigVal->collisionStep && !pConfigVal->output.originDumpPath.empty()) {
        // output atoms in system before collision.
        output(step, true); // dump atoms
    }
}

void MDSimulation::postStep(const unsigned long step) {
    // output atoms information if it is dumping step
    if ((step + 1) % pConfigVal->output.atomsDumpInterval == 0) {
        output(step + 1);
    }
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

#endif //MD_DEV_MODE

void MDSimulation::output(size_t time_step, bool before_collision) {
    // atom boundary in array.
    _type_lattice_coord begin[DIMENSION] = {
            _p_domain->dbx_lattice_coord_sub_box_region.x_low - _p_domain->dbx_lattice_coord_ghost_region.x_low,
            _p_domain->dbx_lattice_coord_sub_box_region.y_low - _p_domain->dbx_lattice_coord_ghost_region.y_low,
            _p_domain->dbx_lattice_coord_sub_box_region.z_low - _p_domain->dbx_lattice_coord_ghost_region.z_low};
    _type_lattice_coord end[DIMENSION] = {
            begin[0] + _p_domain->dbx_lattice_size_sub_box[0],
            begin[1] + _p_domain->dbx_lattice_size_sub_box[1],
            begin[2] + _p_domain->dbx_lattice_size_sub_box[2]};
    _type_lattice_size atoms_size = _p_domain->dbx_lattice_size_sub_box[0] * _p_domain->dbx_lattice_size_sub_box[1] *
                                    _p_domain->dbx_lattice_size_sub_box[2];
    double start = 0, stop = 0;
    static double totalDumpTime = 0;

    start = MPI_Wtime();
    if (!pConfigVal->output.outByFrame) {
        static AtomDump *dumpInstance = nullptr; // pointer used for non-by-frame dumping.
        if (dumpInstance == nullptr) { // initialize atomDump if it is not initialized.
            dumpInstance = new AtomDump(pConfigVal->output.atomsDumpMode, pConfigVal->output.atomsDumpFilePath,
                                        begin, end, atoms_size); // atoms dump.
            // fixme Attempting to use an MPI routine after finalizing MPICH.
        }
        dumpInstance->dump(_atom->getAtomList(), _atom->getInterList(), time_step);
        if (time_step + pConfigVal->output.atomsDumpInterval > pConfigVal->timeSteps) { // the last time of dumping.
            dumpInstance->writeDumpHeader();
            delete dumpInstance;
        }
    } else {
        std::string filename;
        if (before_collision) {
            filename = pConfigVal->output.originDumpPath; // todo pass file name from func output parameters.
        } else {
            filename = fmt::format(pConfigVal->output.atomsDumpFilePath, time_step);
        }
        // pointer to the atom dump class for outputting atoms information.
        auto *dumpInstance = new AtomDump(pConfigVal->output.atomsDumpMode, filename,
                                          begin, end, atoms_size);
        dumpInstance->dump(_atom->getAtomList(), _atom->getInterList(), time_step);
        dumpInstance->writeDumpHeader();
        delete dumpInstance;
    }
    stop = MPI_Wtime();
    totalDumpTime += (stop - start);

    // log dumping time.
    if (time_step + pConfigVal->output.atomsDumpInterval > pConfigVal->timeSteps &&
        MPIDomain::sim_processor.own_rank == MASTER_PROCESSOR) {
        if (pConfigVal->output.atomsDumpMode == OUTPUT_COPY_MODE) {
            kiwi::logs::i("dump", "time of dumping atoms in copy mode:{}.\n", totalDumpTime);
        } else if (pConfigVal->output.atomsDumpMode == OUTPUT_DIRECT_MODE) {
            kiwi::logs::i("dump", "time of dumping atoms in direct mode:{}.\n", totalDumpTime);
        }
    }
}
