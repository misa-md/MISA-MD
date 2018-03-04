#ifndef INPUT_H_
#define INPUT_H_

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

#endif /*INPUT_H_*/
