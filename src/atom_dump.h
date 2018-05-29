//
// Created by genshen on 5/6/18.
//

#ifndef CRYSTAL_MD_ATOM_DUMP_H
#define CRYSTAL_MD_ATOM_DUMP_H

#include "atom.h"
#include "config_values.h"

/**
 * Dump atoms information(including position and velocity of atoms) to binary or text file(s).
 */
class AtomDump {
public:
    AtomDump();

    AtomDump(_type_out_mode mode, const std::string &filename,
             _type_lattice_coord begin[DIMENSION], _type_lattice_coord end[DIMENSION],
             _type_lattice_size atoms_size);

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
    void dump(atom *atom);

private:
    std::string dump_file_name;
    _type_out_mode _dump_mode;

    _type_lattice_size _atoms_size;
    _type_lattice_coord _begin[DIMENSION];
    _type_lattice_coord _end[DIMENSION];

    kiwi::IOWriter *dump_writer = nullptr; // io writer for writing a shared file using mpi-IO lib.

    /**
     * dump atoms with copy mode.
     */
    void dumpModeCopy(atom *atom);

    /*
     * dump atoms with direct mode.
     */
    void dumpModeDirect(atom *atom);
};


#endif //CRYSTAL_MD_ATOM_DUMP_H
