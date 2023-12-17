//
// Created by genshen on 5/6/18.
//

#ifndef MISA_MD_ATOM_DUMP_H
#define MISA_MD_ATOM_DUMP_H

#include <io/local_storage.h>
#include "atom.h"
#include "../config_values.h"
#include "buffered_io.h"
#include "plugin_api.h"

/**
 * Dump atoms information(including position and velocity of atoms) to binary or text file(s).
 */
class AtomDump {
public:
    AtomDump();

    /**
     * dumping atoms information to binary file, including intel atoms.
     * @param frames total frames to be dumped.
     * @param begin starting atoms index of dumping in 3d
     * @param end ending atoms index of dumping in 3d
     * @param atoms_size the count of atoms to be dumped in one frame.
     */
    AtomDump(_type_lattice_size atoms_size, unsigned int frames, atom_dump::type_dump_mask dump_mask,
             _type_lattice_coord begin[DIMENSION], _type_lattice_coord end[DIMENSION]);

    ~AtomDump();

    /**
     * Try to create local storage for storing atoms data.
     * if the local storage is already created, then skip creation.
     * @param dump_file_name the path of dumping file.
     */
    void tryCreateLocalStorage(std::string dump_file_name);

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
     * set header for one frame
     * it should be called before `dumpFrame`.
     * @param time_step current time step.
     */
    void setFrameHeader(size_t time_step);

    /**
     * check IO plugins and dump region.
     * @param region_enabled enabled/disable region dump.
     * @param region the dump region.
     * @param x position of the atom for dumping.
     * @param io_plugins the io plugin for atoms filter.
     * @return true for can dump, false for not dump.
     */
    static bool can_atom_dump(const bool region_enabled, const comm::Region<_type_atom_location> region,
                       const _type_atom_location x[DIMENSION], plugins::IOPlugin *io_plugins);

    /**
     * dumpFrame dump one frame of atoms in current system into file(s), as well as frame header.
     * @param region the dump region
     * @param region_enabled enabled region dump.
     * @param io_plugins the io plugin for atoms filter.
     */
    void dumpFrame(const comm::Region<double> region, const bool region_enabled,
                   AtomList *atom_list, InterAtomList *inter_list, size_t time_step, plugins::IOPlugin *io_plugins);

    /**
     * It write global header and frame headers into dump file.
     */
    void writeDumpHeader();

    /**
     * and close the shared file.
     */
    void onClose();

private:
    const unsigned int frames; // total frame
    unsigned int cur_frame; // current frame
    _type_lattice_size _atoms_size;
    _type_lattice_coord _begin[DIMENSION];
    _type_lattice_coord _end[DIMENSION];

    const atom_dump::type_dump_mask mask;

    std::vector<atom_dump::FrameMetaData> frames_meta;

    kiwi::LocalStorage *local_storage = nullptr; // io writer for writing a shared file using mpi-IO lib.
    BufferedFileWriter *buffered_writer = nullptr;
    MPI_File pFile = NULL; // used in copy mode. // todo close file.
};


#endif //MISA_MD_ATOM_DUMP_H
