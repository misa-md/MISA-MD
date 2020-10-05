//
// Created by genshen on 5/6/18.
//

#ifndef MISA_MD_ATOM_DUMP_H
#define MISA_MD_ATOM_DUMP_H

#include <io/local_storage.h>
#include "atom.h"
#include "../config_values.h"
#include "buffered_io.h"

/**
 * Dump atoms information(including position and velocity of atoms) to binary or text file(s).
 */
class AtomDump {
public:
    AtomDump();

    /**
     * dumping atoms information to binary file, including intel atoms.
     * @param mode dumping mode, mode or direct.
     * @param filename the path of dumping file.
     * @param begin starting atoms index of dumping in 3d
     * @param end ending atoms index of dumping in 3d
     * @param atoms_size the count of atoms to be dumped
     */
    AtomDump(const std::string &filename, _type_lattice_size atoms_size,
             _type_lattice_coord begin[DIMENSION], _type_lattice_coord end[DIMENSION]);

    ~AtomDump();

    /**
     *
     * @param begin  // todo document.
     * @param end
     * @param atoms_size the count of atom to be output.
     * @return
     */
    AtomDump &setBoundary(_type_lattice_coord begin[DIMENSION], _type_lattice_coord end[DIMENSION],
                          _type_lattice_size atoms_size);

    /**
     * set dump file name to store information of dumped atoms.
     * @param filename
     * @return
     */
    AtomDump &setDumpFile(const std::string &filename);

    /**
     * dump atoms to file(s).
     */
    void dump(AtomList *atom_list, InterAtomList *inter_list, size_t time_step);

    void writeDumpHeader();

private:
    std::string _dump_file_name;

    _type_lattice_size _atoms_size;
    _type_lattice_coord _begin[DIMENSION];
    _type_lattice_coord _end[DIMENSION];

    kiwi::LocalStorage *local_storage = nullptr; // io writer for writing a shared file using mpi-IO lib.
    BufferedFileWriter *buffered_writer = nullptr;
    MPI_File pFile = NULL; // used in copy mode. // todo close file.
};


#endif //MISA_MD_ATOM_DUMP_H
