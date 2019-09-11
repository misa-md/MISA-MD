//
// Created by genshen on 2019-07-29.
//

#include <utils/mpi_domain.h>
#include "output_copy.h"

OutputCopy::OutputCopy(const Output output, const comm::BccDomain p_domain)
        : OutputBaseInterface(output, p_domain) {}

void OutputCopy::prepareOutput(const comm::BccDomain domain) {
    // atom boundary in array.
    atoms_size = domain.dbx_sub_box_lattice_size[0] *
                 domain.dbx_sub_box_lattice_size[1] *
                 domain.dbx_sub_box_lattice_size[2];
}

void OutputCopy::onOutputStep(const unsigned long time_step, AtomList *atom_list,
                              InterAtomList *inter_atom_list) {
    double start = 0, stop = 0;
    start = MPI_Wtime();
    if (!output_config.outByFrame) {
        if (dumpInstance == nullptr) { // initialize atomDump if it is not initialized.
            dumpInstance = new AtomDump(output_config.atomsDumpFilePath, atoms_size, begin, end); // atoms dump.
            // fixme Attempting to use an MPI routine after finalizing MPICH.
        }
        dumpInstance->dump(atom_list, inter_atom_list, time_step);
    } else {
        std::string filename = fmt::format(output_config.atomsDumpFilePath, time_step);
        // pointer to the atom dump class for outputting atoms information.
        auto *dumpInstance = new AtomDump(filename, atoms_size, begin, end);
        dumpInstance->dump(atom_list, inter_atom_list, time_step);
        dumpInstance->writeDumpHeader();
        delete dumpInstance;
    }
    stop = MPI_Wtime();
    totalDumpTime += (stop - start);
}

void OutputCopy::beforeCollision(const unsigned long time_step, AtomList *atom_list,
                                 InterAtomList *inter_atom_list) {
    std::string filename = output_config.originDumpPath; // todo pass file name from func output parameters.
    // pointer to the atom dump class for outputting atoms information.
    auto *dumpInstance = new AtomDump(filename, atoms_size, begin, end);
    dumpInstance->dump(atom_list, inter_atom_list, time_step);
    dumpInstance->writeDumpHeader();
    delete dumpInstance;
}

void OutputCopy::onAllOut(const unsigned long time_step) {
    if (!output_config.outByFrame) {
        // if it is not outputting by frame, write header after the final step.
        dumpInstance->writeDumpHeader();
        delete dumpInstance;
    }
    // log dumping time.
    kiwi::logs::i(MASTER_PROCESSOR, "dump", "time of dumping atoms in copy mode: {}.\n", totalDumpTime);
}
