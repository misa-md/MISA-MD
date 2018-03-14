//
// Created by genshen(genshenchu@gmail.com) on 2018-3-13.
//

#ifndef CRYSTAL_MD_SYS_DUMP_H
#define CRYSTAL_MD_SYS_DUMP_H

#include <string>
#include "mpi_utils.h"
#include "data_def.h"


/*
 * This file is basically for writing atom information and system information to file, a shared file.
 * Each processor handle some blocks, which are located by processors rank_id order.
 * Each block have DEFAULT_IO_BLOCK_SIZE bytes. it has a head with DEFAULT_IO_HEADER_SIZE byte data in that file.
 * example: space divide for 4 processors.
 * ---------------------------------------------------
 *|      |    |    |    |    |    |    |    |    |    |
 *| head | p1 | p2 | p3 | p4 | p1 | p2 | p3 | p4 | ...|
 *|      |    |    |    |    |    |    |    |    |    |
 * ---------------------------------------------------
 */

#define DEFAULT_IO_HEADER_SIZE (128) // todo
#define DEFAULT_IO_BLOCK_SIZE (1024*1024) // 1MiB

class IOWriter {
public:

    /**
     * initial variable pDumpFile, and set MPI_IO file view.
     */
    IOWriter(const std::string &filename);

    IOWriter(const std::string &filename, long headerSize, long blockSize);

    /**
     * close file, release unnecessary variable.
     */
    ~IOWriter();

    /**
      *  write data indicated by b to file pDumpFile.
      * @param b  data to be writen.
      * @param start start position of array b.
      * @param size length of array b.
      * @return size that has been writen.
      */// start from 0.
    template<typename T>
    long write(T *b, long start, long size);

    // write from index 0 of array b.
    template<typename T>
    long write(T *b, long size);

protected:
    MPI_File pFile;
private:
    const long blockSize;
    const long headerSize;
};


// methods implements.
template<typename T>
long IOWriter::write(T *b, long start, long size) {
    MPI_File_write(pFile, b + start, size * sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);
    return size;
}

template<typename T>
long IOWriter::write(T *b, long size) {
    MPI_File_write(pFile, b, size * sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);
    return size;
}


#endif //CRYSTAL_MD_SYS_DUMP_H
