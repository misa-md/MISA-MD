//
// Created by baihe back to 2015-08-17.
//

#ifndef MISA_MD_INPUT_H
#define MISA_MD_INPUT_H

#include <string>
#include <fstream>

#include "atom.h"

struct inp_header {
    unsigned long uniq;
    unsigned int version;
    long atom_count;
    long mask;
};

class input {
public:
    input();

    ~input();

    void readPhaseSpace(const std::string read_inp_path, atom *_atom, comm::BccDomain *p_domain);

private:
    static inp_header readHeader(std::fstream &fs);

    /**
     * read atoms from @param fs.
     * @param fs binary file stream.
     * @param head head part of the binary head.
     * @param _atom atom container to store the read atoms.
     * @param p_domain simulation domain.
     * @return the number of atoms read by current process (only count atoms belonging to current process).
     */
    static _type_atom_count readAtoms(std::fstream &fs, const inp_header head, atom *_atom, comm::BccDomain *p_domain);

    static void checkAtomRead(const _type_atom_count atom_added, comm::BccDomain *p_domain);
};

#endif //MISA_MD_INPUT_H
