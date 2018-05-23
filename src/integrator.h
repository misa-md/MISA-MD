//
// Created by baihe back to 2015-05-13.
//

#ifndef CRYSTAL_MD_INTEGRATOR_H
#define CRYSTAL_MD_INTEGRATOR_H

#include "atom.h"

class integrator {
public:
    integrator(double timestepLength);

    ~integrator();

    void firststep(atom *_atom);

    void secondstep(atom *_atom);

    void setTimestepLength(double dt) {
        _timestepLength = dt;
    }

private:
    double _timestepLength;
};

#endif // CRYSTAL_MD_INTEGRATOR_H
