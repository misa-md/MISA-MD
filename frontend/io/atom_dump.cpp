//
// Created by genshen on 2018-05-06.
//

#include <cstdio>
#include <fstream>
#include <logs/logs.h>
#include "atom_dump.h"
#include "utils/mpi_domain.h"
#include "atom_info_dump.h"
#include "types/pre_define.h"

AtomDump::AtomDump() : _dump_file_name(DEFAULT_OUTPUT_DUMP_FILE_PATH), _atoms_size(0) {}

AtomDump::AtomDump(const std::string &filename, _type_lattice_size atoms_size,
                   _type_lattice_coord *begin, _type_lattice_coord *end)
        : _dump_file_name(filename), _atoms_size(atoms_size) {
    setBoundary(begin, end, atoms_size); // todo set variable atoms_size twice.
}

AtomDump::~AtomDump() {
    delete local_storage;
    delete buffered_writer;
    if (pFile != NULL) {
        MPI_File_close(&pFile);
    }
}

AtomDump &AtomDump::setBoundary(_type_lattice_coord *begin, _type_lattice_coord *end, _type_lattice_size atoms_size) {
    for (int i = 0; i < DIMENSION; i++) {
        _begin[i] = begin[i];
        _end[i] = end[i];
    }
    _atoms_size = atoms_size;
    return *this;
}

// todo asynchronous io.
void AtomDump::dump(AtomList *atom_list, InterAtomList *inter_list, size_t time_step) {
    // initialize kiwi writer for copy mode dump.
    if (local_storage == nullptr) {
        int status = MPI_File_open(MPIDomain::sim_processor.comm, _dump_file_name.c_str(), // todo comm.
                                   MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &pFile);
        if (status == -1) {
            kiwi::logs::e("dump", "can not access file {} in MPI IO.", _dump_file_name);
            return; // todo exit application.
        }

        // atom_dump::registerAtomDumpMPIDataType(); // fixme: MPI_Type_contiguous(count=32768, INVALID DATATYPE,
        // make block size equals to buffer size to fix mpi-io writing problems.
        local_storage = new kiwi::LocalStorage(pFile, 128, 128,
                                               BufferedFileWriter::DefaultBufferSize * sizeof(atom_dump::AtomInfoDump));
        local_storage->make(MPI_BYTE, MPIDomain::sim_processor); // create file.
        buffered_writer = new BufferedFileWriter(local_storage, BufferedFileWriter::DefaultBufferSize);
    }

    // dumping inter atoms.
    kiwi::logs::v("dump", "inter atoms count: {}\n", inter_list->nLocalInter());
    for (AtomElement &inter_ref :inter_list->inter_list) {
        buffered_writer->write(&inter_ref, time_step);
    }
    // dumping normal atoms
    for (int k = _begin[2]; k < _end[2]; k++) {
        for (int j = _begin[1]; j < _end[1]; j++) {
            for (int i = _begin[0]; i < _end[0]; i++) {
                AtomElement &atom_ = atom_list->getAtomEleByGhostIndex(i, j, k);
                if (atom_.type == atom_type::INVALID) {
                    continue; // skip out of boxed atoms.
                }
                buffered_writer->write(&atom_, time_step);
            }
        }
    }
    buffered_writer->flush();
}

void AtomDump::writeDumpHeader() {
    struct {
        size_t atoms_count;
    } header{buffered_writer->totalAtomsWritten()};

    local_storage->writeHeader(reinterpret_cast<kiwi::byte *>(&header), sizeof(header), MPIDomain::sim_processor);
}
