//
// Created by genshen on 2019-07-29.
//

#include <fstream>
#include <utils/mpi_domain.h>
#include <logs/logs.h>
#include "output_dump.h"

OutputDump::OutputDump(const DumpConfig output, const comm::BccDomain p_domain)
        : OutputBaseInterface(output, p_domain) {}

void OutputDump::prepareOutput(const comm::BccDomain p_domain) {}

void OutputDump::onOutputStep(const unsigned long time_step, AtomList *atom_list, InterAtomList *inter_atom_list,
                              plugins::IOPlugin *io_plugins) {
  double start = 0, stop = 0;
    start = MPI_Wtime();
    dump(time_step, atom_list, inter_atom_list);
    stop = MPI_Wtime();
    total_dump_time += (stop - start);
}

void OutputDump::onAllOut(const unsigned long time_step) {
    kiwi::logs::i(MASTER_PROCESSOR, "dump", "time of dumping atoms in direct mode: {}.\n", total_dump_time);
}

void OutputDump::dump(const unsigned long time_step, AtomList *atom_list, InterAtomList *inter_atom_list) {
    char outfileName[20];
    sprintf(outfileName, "dump_%d_%ld.atom", MPIDomain::sim_processor.own_rank, time_step);

    std::ofstream outfile;
    outfile.open(outfileName);

    outfile << "print atoms" << std::endl;

    atom_list->foreachSubBoxAtom(
            [atom_list, &outfile](const _type_atom_index gid) {
                MD_LOAD_ATOM_VAR(_atom_ref, atom_list, gid);
                if (!MD_IS_ATOM_TYPE_INTER(_atom_ref, gid)) {
                    outfile << MD_GET_ATOM_ID(_atom_ref, gid) << " "
                            // << "ty" << atom_.type << " "
                            << MD_GET_ATOM_X(_atom_ref, gid, 0) << " "
                            << MD_GET_ATOM_X(_atom_ref, gid, 1) << " "
                            << MD_GET_ATOM_X(_atom_ref, gid, 2) << std::endl;
                }
            }
    );
    outfile << "print inter" << std::endl;
    for (AtomElement &inter_ref :inter_atom_list->inter_list) {
        outfile << inter_ref.id << " "
                // << "ty" << atom->typeinter[i] << " "
                << inter_ref.x[0] << " "
                << inter_ref.x[1] << " "
                << inter_ref.x[2] << std::endl;
    }
    outfile.close();
}
