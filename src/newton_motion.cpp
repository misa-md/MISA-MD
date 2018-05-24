//
// Created by genshen on 5/13/18.
//

#include "newton_motion.h"

NewtonMotion::NewtonMotion(double timestepLength) : _timestepLength(timestepLength) {

}

NewtonMotion::~NewtonMotion() = default;

void NewtonMotion::firststep(AtomList *atom_list, InterAtomList *inter_atom_list) {
    double _m = 55.845;
    double dt_halve = 0.5 * _timestepLength * ftm2v;
    double dtInv2m = dt_halve / _m;  //  1/2 * ftm2v * dt / * _m
    computeFirst(dtInv2m, atom_list, inter_atom_list);
}


void NewtonMotion::secondstep(AtomList *atom_list, InterAtomList *inter_atom_list) {
    double _m = 55.845;
    double dt_halve = 0.5 * _timestepLength * ftm2v;
    double dtInv2m = dt_halve / _m;
    computeSecond(dtInv2m, atom_list, inter_atom_list);
}

void NewtonMotion::computeFirst(double dtInv2m, AtomList *atom_list, InterAtomList *inter_atom_list) {
    double &dt = _timestepLength;
//本地晶格点上的原子求解运动方程第一步
    atom_list->foreachSubBoxAtom(
            [dtInv2m, dt](AtomElement &_atom_ref) {
                if (!_atom_ref.isInterElement()) {
                    _atom_ref.v[0] = _atom_ref.v[0] + dtInv2m * _atom_ref.f[0];
                    _atom_ref.v[1] = _atom_ref.v[1] + dtInv2m * _atom_ref.f[1];
                    _atom_ref.v[2] = _atom_ref.v[2] + dtInv2m * _atom_ref.f[2];
                    _atom_ref.x[0] += dt * _atom_ref.v[0];
                    _atom_ref.x[1] += dt * _atom_ref.v[1];
                    _atom_ref.x[2] += dt * _atom_ref.v[2];
                }
            }
    );

//本地间隙原子求解运动方程第一步
    for (int i = 0; i < inter_atom_list->nlocalinter; i++) {
        for (unsigned short d = 0; d < 3; ++d) {
            inter_atom_list->vinter[i][d] = inter_atom_list->vinter[i][d] + dtInv2m * inter_atom_list->finter[i][d];
            inter_atom_list->xinter[i][d] += dt * inter_atom_list->vinter[i][d];
        }
    }
}


void NewtonMotion::computeSecond(double dtInv2m, AtomList *atom_list, InterAtomList *inter_atom_list) {
    //本地晶格点上的原子求解运动方程第二步
    atom_list->foreachSubBoxAtom(
            [dtInv2m](AtomElement &_atom_ref) {
                // fixme, excluding InterElement. (add if: isInterElement)
                for (unsigned short d = 0; d < DIMENSION; ++d) {
                    _atom_ref.v[d] += dtInv2m * _atom_ref.f[d];
                }
            }
    );
    //本地间隙原子求解运动方程第二步
    for (int i = 0; i < inter_atom_list->nlocalinter; i++) {
        for (unsigned short d = 0; d < 3; ++d) {
            inter_atom_list->vinter[i][d] = inter_atom_list->vinter[i][d] + dtInv2m * inter_atom_list->finter[i][d];
        }
    }
}
