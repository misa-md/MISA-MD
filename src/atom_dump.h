//
// Created by genshen on 5/6/18.
//

#ifndef CRYSTAL_MD_ATOM_DUMP_H
#define CRYSTAL_MD_ATOM_DUMP_H

#include <io/local_storage.h>
#include "atom.h"
#include "config_values.h"

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
    AtomDump(_type_out_mode mode, const std::string &filename,
             _type_lattice_coord begin[DIMENSION], _type_lattice_coord end[DIMENSION],
             _type_lattice_size atoms_size);

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
     * set dump mode, copy or direct(default).
     * @param mode  OUTPUT_DIRECT_MODE or OUTPUT_COPY_MODE
     * @return
     */
    AtomDump &setMode(_type_out_mode mode);

    /**
     * set dump file name to store information of dumped atoms.
     * @param filename
     * @return
     */
    AtomDump &setDumpFile(const std::string &filename);

    /**
     * dump atoms to file(s).
     */
    void dump(atom *atom, size_t time_step);

    void writeDumpHeader();

    void dumpInterLists(InterAtomList *pList);

private:
    std::string _dump_file_name;
    _type_out_mode _dump_mode;

    size_t atom_total = 0; // the count of atoms have writen to local storage.

    _type_lattice_size _atoms_size;
    _type_lattice_coord _begin[DIMENSION];
    _type_lattice_coord _end[DIMENSION];

    kiwi::LocalStorage *local_storage = nullptr; // io writer for writing a shared file using mpi-IO lib.
    MPI_File pFile; // used in copy mode.

    /**
     * dump atoms with copy mode.
     */
    void dumpModeCopy(atom *atom, size_t time_step);

    /*
     * dump atoms with direct mode.
     */
    void dumpModeDirect(atom *atom, size_t time_step);
};


#endif //CRYSTAL_MD_ATOM_DUMP_H
