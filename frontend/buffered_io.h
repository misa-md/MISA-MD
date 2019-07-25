//
// Created by genshen on 2019-07-25.
//

#ifndef CRYSTAL_MD_BUFFERED_FILE_WRITER_H
#define CRYSTAL_MD_BUFFERED_FILE_WRITER_H


#include <mpi.h>
#include <io/local_storage.h>
#include "atom/atom_element.h"
#include "atom_info_dump.h"

class BufferedFileWriter {
public:
    static const unsigned long DefaultBufferSize = 1024;

    /**
     * create a buffered file writer using MPI-IO.
     * In the constructor, we allocate buffer array and
     * set local storage pointer for writing shared file.
     *
     * @param p_file MPI-IO file pointer.
     * @param buffer_size buffer size.
     */
    BufferedFileWriter(kiwi::LocalStorage *p_local_storage, const unsigned long buffer_size);

    ~BufferedFileWriter();

    /**
     * Add atoms to buffer.
     * If the buffer is full, we write buffered data to MPI-IO file and then write left atoms to buffer.
     * @param atom atom pointer.
     * @param time_step current simulation step.
     */
    void write(AtomElement *atom, const unsigned long time_step);

    /**
     * flush atoms data in buffer to MPI-IO file, and clear buffer array.
     */
    void flush();

    /**
     * @return the total atoms have been written.
     */
    inline unsigned long totalAtomsWritten() {
        return total_atoms_written;
    }

private:
    /**
     * io writer for writing a shared file using mpi-IO lib.
     */
    kiwi::LocalStorage *p_local_storage;

    /**
     * array to save buffered atom data.
     */
    atom_dump::AtomInfoDump *buffer;

    /**
     * the size of buffer array.
     */
    const unsigned long buffer_size;

    /**
     * the atoms data have been added into buffer array.
     * (or index of buffer array next atom data will store)
     */
    unsigned long next_buffer_index = 0;

    /**
     * the total atoms have been written.
     */
    unsigned long total_atoms_written = 0;
};


#endif //CRYSTAL_MD_BUFFERED_FILE_WRITER_H
