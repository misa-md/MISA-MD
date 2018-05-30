//
// Created by genshen on 2018-05-06.
//

#include <cstdio>
#include <fstream>
#include <logs/logs.h>
#include "atom_dump.h"

AtomDump::AtomDump() : dump_file_name(DEFAULT_OUTPUT_DUMP_FILE_PATH),
                       _dump_mode(OUTPUT_DIRECT_MODE), _atoms_size(0) {}

AtomDump::AtomDump(_type_out_mode mode, const std::string &filename, _type_lattice_coord *begin,
                   _type_lattice_coord *end, _type_lattice_size atoms_size)
        : _dump_mode(mode), _atoms_size(atoms_size), dump_file_name(filename) {
    setBoundary(begin, end, atoms_size); // todo set variable atoms_size twice.
}

AtomDump::~AtomDump() {
    delete local_storage;
}

AtomDump &AtomDump::setMode(_type_out_mode mode) {
    this->_dump_mode = mode;
    return *this;
}

AtomDump &AtomDump::setDumpFile(const std::string &filename) {
    this->dump_file_name = filename;
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

void AtomDump::dump(atom *atom, size_t time_step) {
    double start, stop;
    if (_dump_mode == OUTPUT_COPY_MODE) { // todo copy atoms, then write.
        start = MPI_Wtime();
        dumpModeCopy(atom, time_step);
        stop = MPI_Wtime();
        kiwi::logs::i("dump", "time of dumping atoms in copy mode:{}.\n", stop - start);
    } else {
        start = MPI_Wtime();
        dumpModeDirect(atom, time_step);
        stop = MPI_Wtime();
        kiwi::logs::i("dump", "time of dumping atoms in direct mode:{}.\n", stop - start);
    }
}

// todo asynchronous io.
void AtomDump::dumpModeCopy(atom *atom, size_t time_step) {
    // initialize kiwi writer for copy mode dump.
    if (local_storage == nullptr) {
        local_storage = new kiwi::LocalStorage(dump_file_name);
        local_storage->make(); // create file.
    }
    long kk;

    const size_t list_buffer_size = 4096;
    size_t list_buff_index = 0;
    struct {
        _type_atom_id id;
        size_t step;
        _type_atom_type type;
//        bool is_inter;
        _type_atom_location atom_location[DIMENSION]; // atom location
        _type_atom_velocity atom_velocity[DIMENSION]; // atom velocity
    } atom_list_buffer[list_buffer_size]; // 4k

    for (int k = _begin[2]; k < _end[2]; k++) {
        for (int j = _begin[1]; j < _end[1]; j++) {
            for (int i = _begin[0]; i < _end[0]; i++) {
                kk = atom->atom_list->IndexOf3DIndex(i, j, k);
                AtomElement &atom_ = atom->getAtomList()->getAtomEleByLinearIndex(kk);
                atom_list_buffer[list_buff_index].id = atom_.id;
                atom_list_buffer[list_buff_index].step = time_step;
                atom_list_buffer[list_buff_index].type = atom_.type;
                atom_list_buffer[list_buff_index].atom_location[0] = atom_.x[0];
                atom_list_buffer[list_buff_index].atom_location[1] = atom_.x[1];
                atom_list_buffer[list_buff_index].atom_location[2] = atom_.x[2];
                atom_list_buffer[list_buff_index].atom_velocity[0] = atom_.v[0];
                atom_list_buffer[list_buff_index].atom_velocity[1] = atom_.v[1];
                atom_list_buffer[list_buff_index].atom_velocity[2] = atom_.v[2];
                list_buff_index++;
                if (list_buff_index == list_buffer_size) { // write data and reset
                    list_buff_index = 0;
                    atom_total += list_buffer_size;
                    local_storage->write(atom_list_buffer, list_buffer_size);
                }
            }
        }
    }
    if (list_buff_index != 0) {
        atom_total += list_buff_index;
        local_storage->write(atom_list_buffer, list_buff_index); // write left data.
    }
}

void AtomDump::dumpModeDirect(atom *atom, size_t time_step) {
    char outfileName[20];
    sprintf(outfileName, "dump_%d_%ld.atom", kiwi::mpiUtils::own_rank, time_step);

    std::ofstream outfile;
    outfile.open(outfileName);

    outfile << "print atoms" << std::endl;

    atom->atom_list->foreachSubBoxAtom(
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
    for (int i = 0; i < atom->inter_atom_list->nlocalinter; i++) {
        outfile << atom->inter_atom_list->idinter[i] << " "
                //                << "ty" << atom->typeinter[i] << " "
                << atom->inter_atom_list->xinter[i][0] << " "
                << atom->inter_atom_list->xinter[i][1] << " "
                << atom->inter_atom_list->xinter[i][2] << std::endl;
    }
    outfile.close();
}

void AtomDump::writeDumpHeader() {
    struct {
        size_t atoms_count;
    } header{atom_total};

    local_storage->writeHeader(reinterpret_cast<kiwi::byte *>(&header), sizeof(header));
}
