//
// Created by genshen on 5/6/18.
//

#include <cstdio>
#include <fstream>
#include <logs/logs.h>
#include "atom_dump.h"

AtomDump::AtomDump() : dump_file_name(DEFAULT_OUTPUT_DUMP_FILENAME),
                       _dump_mode(OUTPUT_DIRECT_MODE), _atoms_size(0) {}

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

void AtomDump::dump(atom *atom) {
    double start, stop;
    if (_dump_mode == OUTPUT_COPY_MODE) { // todo copy atoms, then write.
        start = MPI_Wtime();
        dumpModeCopy(atom);
        stop = MPI_Wtime();
        kiwi::logs::i("dump", "time of dumping atoms in copy mode:{}.\n", stop - start);
    } else {
        start = MPI_Wtime();
        dumpModeDirect(atom);
        stop = MPI_Wtime();
        kiwi::logs::i("dump", "time of dumping atoms in direct mode:{}.\n", stop - start);
    }
}

// todo asynchronous io.
void AtomDump::dumpModeCopy(atom *atom) {
    // initialize kiwi writer for copy mode dump.
    if (dump_writer == nullptr) {
        dump_writer = new kiwi::IOWriter(dump_file_name);
    }
    long kk;
    double *x_io = new double[_atoms_size * 4];
//        int fd, ret;
//        fd = open(outfileName, O_CREAT | O_TRUNC | O_RDWR, 0700);
//        if (fd == -1) {
//            printf("ERROR,open file %s failed\n", outfileName);
//            abort(1);
//        }

    int n = 0;
    for (int k = _begin[2]; k < _end[2]; k++) {
        for (int j = _begin[1]; j < _end[1]; j++) {
            for (int i = _begin[0]; i < _end[0]; i++) {
                kk = atom->IndexOf3DIndex(i, j, k);
                AtomElement &atom_ = atom->getAtomList()->getAtomEleByLinearIndex(kk);
                x_io[n * 4] = atom_.id;
                x_io[n * 4 + 1] = atom_.x[0];
                x_io[n * 4 + 2] = atom_.x[1];
                x_io[n * 4 + 3] = atom_.x[2];
                n++;
            }
        }
    }
    dump_writer->write(x_io, _atoms_size * 4);
    delete[] x_io;
}

void AtomDump::dumpModeDirect(atom *atom) {
    char outfileName[20];
    sprintf(outfileName, "dump_%d.atom", kiwi::mpiUtils::own_rank);

    ofstream outfile;
    outfile.open(outfileName);

    outfile << "print atoms" << std::endl;

    long kk = 0;
    for (int k = _begin[2]; k < _end[2]; k++) { // todo int type of k.
        for (int j = _begin[1]; j < _end[1]; j++) {
            for (int i = _begin[0]; i < _end[0]; i++) {
                kk = atom->IndexOf3DIndex(i, j, k);
                AtomElement &atom_ = atom->getAtomList()->getAtomEleByLinearIndex(kk);
                if (!atom_.isInterElement())
                    outfile << atom_.id << " "
//                            << "ty" << atom_.type << " "
                            << atom_.x[0] << " "
                            << atom_.x[1] << " "
                            << atom_.x[2] << std::endl;
            }
        }
    }
    outfile << "print inter" << std::endl;
    for (int i = 0; i < atom->nlocalinter; i++) {
        outfile << atom->idinter[i] << " "
//                << "ty" << atom->typeinter[i] << " "
                << atom->xinter[i][0] << " "
                << atom->xinter[i][1] << " "
                << atom->xinter[i][2] << std::endl;
    }
    outfile.close();
}
