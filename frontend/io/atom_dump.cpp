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

AtomDump::AtomDump() : _dump_file_name(DEFAULT_OUTPUT_DUMP_FILE_PATH),
                       _dump_mode(OUTPUT_DIRECT_MODE), _atoms_size(0) {}

AtomDump::AtomDump(_type_out_mode mode, const std::string &filename, _type_lattice_coord *begin,
                   _type_lattice_coord *end, _type_lattice_size atoms_size)
        : _dump_file_name(filename), _dump_mode(mode), _atoms_size(atoms_size) {
    setBoundary(begin, end, atoms_size); // todo set variable atoms_size twice.
}

AtomDump::~AtomDump() {
    delete local_storage;
    delete buffered_writer;
    if (pFile != NULL) {
        MPI_File_close(&pFile);
    }
}

AtomDump &AtomDump::setMode(_type_out_mode mode) {
    this->_dump_mode = mode;
    return *this;
}

AtomDump &AtomDump::setDumpFile(const std::string &filename) {
    this->_dump_file_name = filename;
    return *this;
}

AtomDump &AtomDump::setBoundary(_type_lattice_coord *begin, _type_lattice_coord *end, _type_lattice_size atoms_size) {
    for (int i = 0; i < DIMENSION; i++) {
        _begin[i] = begin[i];
        _end[i] = end[i];
    }
    _atoms_size = atoms_size;
    return *this;
}

void AtomDump::dump(AtomList *atom_list, InterAtomList *inter_list, size_t time_step) {
    if (_dump_mode == OUTPUT_COPY_MODE) { // todo copy atoms, then write.
        dumpModeCopy(atom_list, inter_list, time_step);
    } else {
        dumpModeDirect(atom_list, inter_list, time_step);
    }
}

// todo asynchronous io.
void AtomDump::dumpModeCopy(AtomList *atom_list, InterAtomList *inter_list, size_t time_step) {
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
    _type_atom_index kk;
    for (int k = _begin[2]; k < _end[2]; k++) {
        for (int j = _begin[1]; j < _end[1]; j++) {
            for (int i = _begin[0]; i < _end[0]; i++) {
                kk = atom_list->lattice.IndexOf3DIndex(i, j, k);
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                if (atom_.type == atom_type::INVALID) {
                    continue; // skip out of boxed atoms.
                }
                buffered_writer->write(&atom_, time_step);
            }
        }
    }
    buffered_writer->flush();
}

void AtomDump::dumpModeDirect(AtomList *atom_list, InterAtomList *inter_list, size_t time_step) {
    char outfileName[20];
    sprintf(outfileName, "dump_%d_%ld.atom", MPIDomain::sim_processor.own_rank, time_step);

    std::ofstream outfile;
    outfile.open(outfileName);

    outfile << "print atoms" << std::endl;

    atom_list->foreachSubBoxAtom(
            [&outfile](AtomElement &_atom_ref) {
                if (!_atom_ref.isInterElement()) {
                    outfile << _atom_ref.id << " "
                            // << "ty" << atom_.type << " "
                            << _atom_ref.x[0] << " "
                            << _atom_ref.x[1] << " "
                            << _atom_ref.x[2] << std::endl;
                }
            }
    );
    outfile << "print inter" << std::endl;
    for (AtomElement &inter_ref :inter_list->inter_list) {
        outfile << inter_ref.id << " "
                // << "ty" << atom->typeinter[i] << " "
                << inter_ref.x[0] << " "
                << inter_ref.x[1] << " "
                << inter_ref.x[2] << std::endl;
    }
    for (int i = 0; i < inter_list->nLocalInter(); i++) {

    }
    outfile.close();
}

void AtomDump::writeDumpHeader() {
    struct {
        size_t atoms_count;
    } header{buffered_writer->totalAtomsWritten()};

    local_storage->writeHeader(reinterpret_cast<kiwi::byte *>(&header), sizeof(header), MPIDomain::sim_processor);
}
