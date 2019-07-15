//
// Created by baihe back to 2015-08-17.
//

#ifndef CRYSTAL_MD_INPUT_H
#define CRYSTAL_MD_INPUT_H

#include <string>
#include <fstream>

#include "atom.h"

class input {
public:
    input();

    ~input();

    void readPhaseSpace(atom *_atom, comm::BccDomain *p_domain);

private:
    std::string _phaseSpaceFile;
    std::string _phaseSpaceHeaderFile;
    std::fstream _phaseSpaceFileStream;
    std::fstream _phaseSpaceHeaderFileStream;
};

#endif //CRYSTAL_MD_INPUT_H
