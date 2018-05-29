//
// Created by genshen on 5/13/18.
//

#include "newton_motion.h"

NewtonMotion::NewtonMotion(double timestepLength) : _timestepLength(timestepLength) {
    preComputeDtInv2m();
}

NewtonMotion::~NewtonMotion() = default;

void NewtonMotion::setTimestepLength(double dt) {
    _timestepLength = dt;
    // refresh dtInv2m.
    preComputeDtInv2m();
}

void NewtonMotion::preComputeDtInv2m() {
    atom_type::atom_type type;
    double _m, dt_halve;
    for (int i = 0; i < atom_type::num_atom_types; i++) {
        type = atom_type::getAtomTypeByNum(i);
        _m = atom_type::getAtomMass(type);
        dt_halve = 0.5 * _timestepLength * ftm2v;
        dt_inv_m[i] = dt_halve / _m;  //  1/2 * ftm2v * dt / * _m
    }
}

void NewtonMotion::firststep(AtomList *atom_list, InterAtomList *inter_atom_list) {

    computeFirst(atom_list, inter_atom_list);
}

void NewtonMotion::secondstep(AtomList *atom_list, InterAtomList *inter_atom_list) {
    computeSecond(atom_list, inter_atom_list);
}

void NewtonMotion::computeFirst(AtomList *atom_list, InterAtomList *inter_atom_list) {
    double &dt = _timestepLength;
    //本地晶格点上的原子求解运动方程第一步
    atom_list->foreachSubBoxAtom(
            [=](AtomElement &_atom_ref) {
                if (!_atom_ref.isInterElement()) {
                    _atom_ref.v[0] = _atom_ref.v[0] + dt_inv_m[_atom_ref.type] * _atom_ref.f[0];
                    _atom_ref.v[1] = _atom_ref.v[1] + dt_inv_m[_atom_ref.type] * _atom_ref.f[1];
                    _atom_ref.v[2] = _atom_ref.v[2] + dt_inv_m[_atom_ref.type] * _atom_ref.f[2];
                    _atom_ref.x[0] += dt * _atom_ref.v[0];
                    _atom_ref.x[1] += dt * _atom_ref.v[1];
                    _atom_ref.x[2] += dt * _atom_ref.v[2];
                }
            }
    );

    //本地间隙原子求解运动方程第一步
    atom_type::atom_type _type;
    for (int i = 0; i < inter_atom_list->nlocalinter; i++) {
        for (unsigned short d = 0; d < 3; ++d) {
            _type = inter_atom_list->typeinter[i];
            inter_atom_list->vinter[i][d] = inter_atom_list->vinter[i][d] + dt_inv_m[_type] * inter_atom_list->finter[i][d];
            inter_atom_list->xinter[i][d] += dt * inter_atom_list->vinter[i][d];
        }
    }
}


void NewtonMotion::computeSecond(AtomList *atom_list, InterAtomList *inter_atom_list) {
    //本地晶格点上的原子求解运动方程第二步
    atom_list->foreachSubBoxAtom(
            [=](AtomElement &_atom_ref) {
                // fixme, excluding InterElement. (add if: isInterElement)
                for (unsigned short d = 0; d < DIMENSION; ++d) {
                    _atom_ref.v[d] += dt_inv_m[_atom_ref.type] * _atom_ref.f[d];
                }
            }
    );
    //本地间隙原子求解运动方程第二步
    atom_type::atom_type _type;
    for (int i = 0; i < inter_atom_list->nlocalinter; i++) {
        for (unsigned short d = 0; d < 3; ++d) {
            _type = inter_atom_list->typeinter[i];
            inter_atom_list->vinter[i][d] =
                    inter_atom_list->vinter[i][d] + dt_inv_m[_type] * inter_atom_list->finter[i][d];
        }
    }
}
