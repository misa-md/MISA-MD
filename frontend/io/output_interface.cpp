//
// Created by genshen on 2019-07-29.
//

#include <utils/mpi_domain.h>
#include "output_interface.h"

OutputInterface::OutputInterface(const Output output) : output_config(output) {}

void OutputInterface::prepareOutput(const comm::BccDomain domain) {
    // atom boundary in array.
    begin[0] = domain.dbx_lattice_coord_sub_box_region.x_low - domain.dbx_lattice_coord_ghost_region.x_low;
    begin[1] = domain.dbx_lattice_coord_sub_box_region.y_low - domain.dbx_lattice_coord_ghost_region.y_low;
    begin[2] = domain.dbx_lattice_coord_sub_box_region.z_low - domain.dbx_lattice_coord_ghost_region.z_low;
    end[0] = begin[0] + domain.dbx_lattice_size_sub_box[0];
    end[1] = begin[1] + domain.dbx_lattice_size_sub_box[1];
    end[2] = begin[2] + domain.dbx_lattice_size_sub_box[2];
    atoms_size = domain.dbx_lattice_size_sub_box[0] *
                 domain.dbx_lattice_size_sub_box[1] *
                 domain.dbx_lattice_size_sub_box[2];
}

void OutputInterface::onOutputStep(const unsigned long time_step, AtomList *atom_list,
                                   InterAtomList *inter_atom_list) {
    double start = 0, stop = 0;
    start = MPI_Wtime();
    if (!output_config.outByFrame) {
        if (dumpInstance == nullptr) { // initialize atomDump if it is not initialized.
            dumpInstance = new AtomDump(output_config.atomsDumpMode, output_config.atomsDumpFilePath,
                                        begin, end, atoms_size); // atoms dump.
            // fixme Attempting to use an MPI routine after finalizing MPICH.
        }
        dumpInstance->dump(atom_list, inter_atom_list, time_step);
    } else {
        std::string filename = fmt::format(output_config.atomsDumpFilePath, time_step);
        // pointer to the atom dump class for outputting atoms information.
        auto *dumpInstance = new AtomDump(output_config.atomsDumpMode, filename,
                                          begin, end, atoms_size);
        dumpInstance->dump(atom_list, inter_atom_list, time_step);
        dumpInstance->writeDumpHeader();
        delete dumpInstance;
    }
    stop = MPI_Wtime();
    totalDumpTime += (stop - start);
}

void OutputInterface::beforeCollision(const unsigned long time_step, AtomList *atom_list,
                                       InterAtomList *inter_atom_list) {
    std::string filename = output_config.originDumpPath; // todo pass file name from func output parameters.
    // pointer to the atom dump class for outputting atoms information.
    auto *dumpInstance = new AtomDump(output_config.atomsDumpMode, filename, begin, end, atoms_size);
    dumpInstance->dump(atom_list, inter_atom_list, time_step);
    dumpInstance->writeDumpHeader();
    delete dumpInstance;
}

void OutputInterface::onAllOut(const unsigned long time_step) {
    if (!output_config.outByFrame) {
        // if it is not outputting by frame, write header after the final step.
        dumpInstance->writeDumpHeader();
        delete dumpInstance;
    }
    // log dumping time.
    if (MPIDomain::sim_processor.own_rank == MASTER_PROCESSOR) {
        if (output_config.atomsDumpMode == OUTPUT_COPY_MODE) {
            kiwi::logs::i("dump", "time of dumping atoms in copy mode:{}.\n", totalDumpTime);
        } else if (output_config.atomsDumpMode == OUTPUT_DIRECT_MODE) {
            kiwi::logs::i("dump", "time of dumping atoms in direct mode:{}.\n", totalDumpTime);
        }
    }
}
