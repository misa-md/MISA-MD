#include "integrator.h"

integrator::integrator(double timestepLength) : _timestepLength(timestepLength) {}

integrator::~integrator() {}

void integrator::firststep(atom *_atom) {
    double _m = 55.845;
    double dt_halve = 0.5 * _timestepLength * ftm2v;
    double dtInv2m = dt_halve / _m;  //  1/2 * ftm2v * dt / * _m
    _atom->computefirst(dtInv2m, _timestepLength);
}

void integrator::secondstep(atom *_atom) {
    double _m = 55.845;
    double dt_halve = 0.5 * _timestepLength * ftm2v;
    double dtInv2m = dt_halve / _m;
    _atom->computesecond(dtInv2m);
}
