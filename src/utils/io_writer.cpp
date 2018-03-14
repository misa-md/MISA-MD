//
// Created by genshen(genshenchu@gmail.com) on 2018-3-13.
//

#include <iostream>
#include "io_writer.h"

IOWriter::IOWriter(const std::string &filename) :
        IOWriter(filename, DEFAULT_IO_HEADER_SIZE, DEFAULT_IO_BLOCK_SIZE) {
}

IOWriter::IOWriter(const std::string &filename, long headerSize, long blockSize) :
        headerSize(headerSize), blockSize(blockSize) {
    int status = MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                               MPI_MODE_CREATE | MPI_MODE_WRONLY,
                               MPI_INFO_NULL, &pFile);
    if (status == -1) { // todo no output.
        std::cout << "ERROR,open file" << filename << "failed" << std::endl;
        exit(1);
    }

    MPI_Aint lb, extent;
    MPI_Datatype etype, contig, filetype;

    etype = MPI_BYTE;
    lb = 0; // DEFAULT_IO_BLOCK_SIZE * mpiUtils::ownRank;
    extent = blockSize * mpiUtils::allRanks;

    MPI_Type_contiguous(blockSize, MPI_BYTE, &contig);
    MPI_Type_create_resized(contig, lb, extent, &filetype);
    MPI_Type_commit(&filetype);
    MPI_File_set_view(pFile, headerSize + blockSize * mpiUtils::ownRank,
                      etype, filetype, "native", MPI_INFO_NULL);
}

IOWriter::~IOWriter() {
    if (pFile != nullptr) {
        MPI_File_close(&pFile);
        pFile = nullptr;
    }
}
