//
// Created by genshen on 2018-05-06.
//

#include <cstdio>
#include <fstream>
#include <logs/logs.h>
#include "atom_dump.h"
#include "utils/mpi_domain.h"
#include "atom_dump_types.h"
#include "types/pre_define.h"

AtomDump::AtomDump() : frames(0), cur_frame(0), _atoms_size(0), frames_meta() {}

AtomDump::AtomDump(_type_lattice_size atoms_size, unsigned int frames,
                   _type_lattice_coord *begin, _type_lattice_coord *end)
        : frames(frames), cur_frame(0), _atoms_size(atoms_size), frames_meta() {
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

void AtomDump::tryCreateLocalStorage(std::string dump_file_name) {

    if (local_storage != nullptr) {
        return;
    }
    if (dump_file_name.empty()) {
        dump_file_name = DEFAULT_OUTPUT_DUMP_FILE_PATH;
    }

    // initialize kiwi writer for copy mode dump.
    int status = MPI_File_open(MPIDomain::sim_processor.comm, dump_file_name.c_str(), // todo comm.
                               MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &pFile);
    if (status == -1) {
        kiwi::logs::e("dump", "can not access file {} in MPI IO.\n", dump_file_name);
        // todo exit application.
        return;
    }

    MPI_Datatype mpi_data_type = atom_dump::registerAtomDumpMPIDataType(mask);
    int data_type_size = 0;
    MPI_Type_size(mpi_data_type, &data_type_size);
    // make block size equals to buffer size to fix mpi-io writing problems.
    local_storage = new kiwi::LocalStorage(pFile, sizeof(atom_dump::GlobalMetaData) +
                                                  frames * sizeof(atom_dump::FrameMetaData),
                                           0, // it will be ignored if we don.t call local_storage.make
                                           data_type_size * BufferedFileWriter::DefaultBufferSize);
    // todo un-commit MPI data type
    // create file view.
    local_storage->make(MPI_BYTE, MPIDomain::sim_processor);
    buffered_writer = new BufferedFileWriter(local_storage, BufferedFileWriter::DefaultBufferSize);
}

void AtomDump::setFrameHeader(size_t time_step) {
    if ((cur_frame % MPIDomain::sim_processor.all_ranks) == MPIDomain::sim_processor.own_rank) {
        this->frames_meta.emplace_back(
                atom_dump::FrameMetaData{.step = time_step});
    }
}

// todo asynchronous io.
void AtomDump::dumpFrame(AtomList *atom_list, InterAtomList *inter_list, size_t time_step) {
    if (local_storage == nullptr) {
        kiwi::logs::e("dump", "local storage is not created");
        // todo exit application.
        return;
    }
    if (cur_frame >= frames) {
        kiwi::logs::e("dump", "frame overflow");
        return;
    }

    const MPI_Datatype mpi_data_type = atom_dump::registerAtomDumpMPIDataType(mask);

    // dumping inter atoms.
    kiwi::logs::v("dump", "inter atoms count: {}\n", inter_list->nLocalInter());
    for (AtomElement &inter_ref :inter_list->inter_list) {
        buffered_writer->write(&inter_ref, mpi_data_type);
    }

    // dumping normal atoms
    for (int k = _begin[2]; k < _end[2]; k++) {
        for (int j = _begin[1]; j < _end[1]; j++) {
            for (int i = _begin[0]; i < _end[0]; i++) {
                AtomElement &atom_ = atom_list->getAtomEleByGhostIndex(i, j, k);
                if (atom_.type == atom_type::INVALID) {
                    continue; // skip out of boxed atoms.
                }
                buffered_writer->write(&atom_, mpi_data_type);
            }
        }
    }

    // write terminator atom
    AtomElement atom_;
    atom_.id = 1024;
    atom_.type = atom_type::atom_type::INVALID;
    buffered_writer->write(&atom_, mpi_data_type);

    buffered_writer->flush(mpi_data_type);
    cur_frame++;
}

void AtomDump::writeDumpHeader() {
    // dump frames header
    const size_t g_header_size = sizeof(atom_dump::GlobalMetaData);
    size_t header_size = this->frames_meta.size();
    local_storage->writer.make(g_header_size, sizeof(atom_dump::FrameMetaData),
                               atom_dump::registerFrameMetaMPIDataType(),
                               MPIDomain::sim_processor);
    local_storage->writer.writeAll(this->frames_meta.data(), header_size);

    // dump global header
    local_storage->writer.make(0, sizeof(atom_dump::GlobalMetaData), MPI_BYTE, MPIDomain::sim_processor);
    if (MPIDomain::sim_processor.own_rank == MASTER_PROCESSOR) {
        MPI_Datatype mpi_data_type = atom_dump::registerAtomDumpMPIDataType(mask);
        int data_type_size = 0;
        MPI_Type_size(mpi_data_type, &data_type_size);

        atom_dump::GlobalMetaData global_meta = atom_dump::GlobalMetaData{
                .self_size =  sizeof(atom_dump::GlobalMetaData),
                .frame_meta_size = sizeof(atom_dump::FrameMetaData),
                .block_atoms = BufferedFileWriter::DefaultBufferSize,
                .atoms_num = 0,
                .atom_item_bytes = static_cast<size_t>(data_type_size),
                .mpi_ranks = static_cast<size_t>(MPIDomain::sim_processor.all_ranks),
                .mask= this->mask,
                .format_version=  1,
                .global_header_size= sizeof(atom_dump::GlobalMetaData),
                .local_size= 0,
                .frames= 1,
        };
        local_storage->writer.write(&global_meta, sizeof(atom_dump::GlobalMetaData));
    }
//    atom_dump::GlobalMetaData header{buffered_writer->totalAtomsWritten(), 0x001, 0};
//    local_storage->writeHeader(reinterpret_cast<kiwi::byte *>(&header), sizeof(header), MPIDomain::sim_processor);
}

void AtomDump::onClose() {
//    local_storage->writer
}
