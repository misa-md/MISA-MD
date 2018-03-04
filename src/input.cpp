#include "input.h"
#include "domaindecomposition.h"
#include <iostream>

using namespace std;

input::input() {}

input::~input() {}

void input::readPhaseSpace(atom *_atom) {
    particledata *sendbuf;
    _phaseSpaceFile = "dump.atom";
    cout << "Opening phase space file " << _phaseSpaceFile << endl;
    _phaseSpaceFileStream.open(_phaseSpaceFile.c_str());
    if (!_phaseSpaceFileStream.is_open()) {
        cout << "Could not open phaseSpaceFile " << _phaseSpaceFile << endl;
    }
    cout << "Reading phase space file " << _phaseSpaceFile << endl;
    double x, y, z, vx, vy, vz;
    unsigned long id;
    int type;
    x = y = z = vx = vy = vz = 0.;
    for (unsigned long i = 0; i < 432; i++) {
        _phaseSpaceFileStream >> id >> x >> y >> z >> vx >> vy >> vz;
        _atom->addatom(id, x, y, z, vx, vy, vz);
    }
    _phaseSpaceFileStream.close();
}
