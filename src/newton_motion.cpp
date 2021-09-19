//
// Created by genshen on 5/13/18.
//

#include "newton_motion.h"

NewtonMotion::NewtonMotion(double timestepLength, const _type_atom_types num_types)
        : num_types(num_types), _timestepLength(timestepLength), dt_inv_m(num_types) {
    preComputeDtInv2m();
}

NewtonMotion::~NewtonMotion() = default;

void NewtonMotion::setTimestepLength(const double dt) {
    _timestepLength = dt;
    // refresh dtInv2m.
    preComputeDtInv2m();
}

void NewtonMotion::preComputeDtInv2m() {
    atom_type::atom_type type;
    double _m, dt_halve;
    for (int i = 0; i < num_types; i++) {
        type = atom_type::getAtomTypeByNum(i);
        _m = atom_type::getAtomMass(type);
        dt_halve = 0.5 * _timestepLength * ftm2v;
        dt_inv_m[i] = dt_halve / _m;  //  1/2 * ftm2v * dt / * _m
    }
}

void NewtonMotion::firststep(AtomList *atom_list, InterAtomList *inter_atom_list) {
    double &dt = _timestepLength;
    //本地晶格点上的原子求解运动方程第一步
    atom_list->foreachSubBoxAtom(
            [=](const _type_atom_index gid) {
                MD_LOAD_ATOM_VAR(_atom_ref, atom_list, gid);
                if (!MD_IS_ATOM_TYPE_INTER(_atom_ref, gid)) {
                    MD_ADD_ATOM_V(_atom_ref, gid, 0,
                                  dt_inv_m[MD_GET_ATOM_TYPE(_atom_ref, gid)] * MD_GET_ATOM_F(_atom_ref, gid, 0));
                    MD_ADD_ATOM_V(_atom_ref, gid, 1,
                                  dt_inv_m[MD_GET_ATOM_TYPE(_atom_ref, gid)] * MD_GET_ATOM_F(_atom_ref, gid, 1));
                    MD_ADD_ATOM_V(_atom_ref, gid, 2,
                                  dt_inv_m[MD_GET_ATOM_TYPE(_atom_ref, gid)] * MD_GET_ATOM_F(_atom_ref, gid, 2));

                    MD_ADD_ATOM_X(_atom_ref, gid, 0, dt * MD_GET_ATOM_V(_atom_ref, gid, 0));
                    MD_ADD_ATOM_X(_atom_ref, gid, 1, dt * MD_GET_ATOM_V(_atom_ref, gid, 1));
                    MD_ADD_ATOM_X(_atom_ref, gid, 2, dt * MD_GET_ATOM_V(_atom_ref, gid, 2));
                }
            }
    );

    //本地间隙原子求解运动方程第一步
    atom_type::atom_type _type;
    for (AtomElement &inter_ref :inter_atom_list->inter_list) {
        for (unsigned short d = 0; d < 3; ++d) {
            _type = inter_ref.type;
            inter_ref.v[d] = inter_ref.v[d] + dt_inv_m[_type] * inter_ref.f[d];
            inter_ref.x[d] += dt * inter_ref.v[d];
        }
    }
}

void NewtonMotion::secondstep(AtomList *atom_list, InterAtomList *inter_atom_list) {
    //本地晶格点上的原子求解运动方程第二步
    atom_list->foreachSubBoxAtom(
            [=](const _type_atom_index gid) {
                MD_LOAD_ATOM_VAR(_atom_ref, atom_list, gid);
                if (!MD_IS_ATOM_TYPE_INTER(_atom_ref, gid)) {
                    MD_ADD_ATOM_V(_atom_ref, gid, 0,
                                  dt_inv_m[MD_GET_ATOM_TYPE(_atom_ref, gid)] * MD_GET_ATOM_F(_atom_ref, gid, 0));
                    MD_ADD_ATOM_V(_atom_ref, gid, 1,
                                  dt_inv_m[MD_GET_ATOM_TYPE(_atom_ref, gid)] * MD_GET_ATOM_F(_atom_ref, gid, 1));
                    MD_ADD_ATOM_V(_atom_ref, gid, 2,
                                  dt_inv_m[MD_GET_ATOM_TYPE(_atom_ref, gid)] * MD_GET_ATOM_F(_atom_ref, gid, 2));
                }
            }
    );
    //本地间隙原子求解运动方程第二步
    for (AtomElement &inter_ref :inter_atom_list->inter_list) {
        inter_ref.v[0] += dt_inv_m[inter_ref.type] * inter_ref.f[0];
        inter_ref.v[1] += dt_inv_m[inter_ref.type] * inter_ref.f[1];
        inter_ref.v[2] += dt_inv_m[inter_ref.type] * inter_ref.f[2];
    }
}
