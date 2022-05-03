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

    static void readAtoms(std::fstream &fs, const inp_header head, atom *_atom, comm::BccDomain *p_domain);
};

#endif //MISA_MD_INPUT_H
