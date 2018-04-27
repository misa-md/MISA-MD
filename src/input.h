//
// Created by baihe back to 2015-08-17.
//

#ifndef CRYSTAL_MD_INPUT_H
#define CRYSTAL_MD_INPUT_H

#include "atom.h"

#include <string>
#include <fstream>

class input {
public:
    input();

    ~input();

    void readPhaseSpace(atom *_atom);

private:
    std::string _phaseSpaceFile;
    std::string _phaseSpaceHeaderFile;
    std::fstream _phaseSpaceFileStream;
    std::fstream _phaseSpaceHeaderFileStream;
};

#endif //CRYSTAL_MD_INPUT_H
