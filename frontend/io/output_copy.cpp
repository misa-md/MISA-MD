//
// Created by genshen on 2019-07-29.
//

#include <utils/mpi_domain.h>
#include <logs/logs.h>
#include "output_copy.h"

OutputCopy::OutputCopy(const DumpConfig output, const comm::BccDomain p_domain)
        : dumpInstance(nullptr), OutputBaseInterface(output, p_domain) {}

void OutputCopy::prepareOutput(const comm::BccDomain domain) {
    // atom boundary in array.
    atoms_size = domain.dbx_sub_box_lattice_size[0] *
                 domain.dbx_sub_box_lattice_size[1] *
                 domain.dbx_sub_box_lattice_size[2];
}

void OutputCopy::onOutputStep(const unsigned long time_step, AtomList *atom_list, InterAtomList *inter_atom_list) {
    double start = 0, stop = 0;
    start = MPI_Wtime();
    if (!output_config.by_frame) {
        if (dumpInstance == nullptr) { // initialize atomDump if it is not initialized.
            dumpInstance = new AtomDump(atoms_size, output_config.steps, output_config.dump_mask, begin,
                                        end); // atoms dump.
            dumpInstance->tryCreateLocalStorage(output_config.file_path);
            // fixme Attempting to use an MPI routine after finalizing MPICH.
        }
        dumpInstance->setFrameHeader(time_step);
        dumpInstance->dumpFrame(region, !(output_config.dump_whole_system),
                                atom_list, inter_atom_list, time_step);
    } else {
        std::string filename = fmt::format(output_config.file_path, time_step);
        // pointer to the atom dump class for outputting atoms information.
        auto *dump_instance = new AtomDump(atoms_size, 1, output_config.dump_mask, begin, end);
        dump_instance->tryCreateLocalStorage(filename);
        dump_instance->setFrameHeader(time_step);
        dump_instance->dumpFrame(region, !(output_config.dump_whole_system), atom_list, inter_atom_list, time_step);
        dump_instance->writeDumpHeader();
        delete dump_instance;
    }
    stop = MPI_Wtime();
    totalDumpTime += (stop - start);
}

void OutputCopy::onAllOut(const unsigned long time_step) {
    if (!output_config.by_frame) {
        // if it is not outputting by frame, write header after the final step.
        dumpInstance->writeDumpHeader();
        delete dumpInstance;
    }
    // log dumping time.
    kiwi::logs::i(MASTER_PROCESSOR, "dump", "time of dumping atoms in copy mode: {}.\n", totalDumpTime);
}
