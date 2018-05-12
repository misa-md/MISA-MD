//
// Created by genshen on 5/5/18.
//

#include <cmath>
#include <iostream>
#include <logs/logs.h>
#include "world_builder.h"

WorldBuilder::WorldBuilder() : _random_seed(1024), tset(0),
                               box_x(0), box_y(0), box_z(0) {
    _p_domain = nullptr;
    _p_atom = nullptr;
}

WorldBuilder &WorldBuilder::setDomain(Domain *p_domain) {
    this->_p_domain = p_domain;
    return *this;
}

WorldBuilder &WorldBuilder::setAtomsContainer(atom *p_atom) {
    this->_p_atom = p_atom;
    return *this;
}

WorldBuilder &WorldBuilder::setRandomSeed(int seed) {
    this->_random_seed = seed;
    return *this;
}

WorldBuilder &WorldBuilder::setTset(double tset) {
    this->tset = tset;
    return *this;
}

WorldBuilder &WorldBuilder::setLatticeConst(double lattice_const) {
    this->_lattice_const = lattice_const;
    return *this;
}

WorldBuilder &WorldBuilder::setAlloyRatio(int ratio[atom_type::num_atom_types]) {
    for (int i = 0; i < atom_type::num_atom_types; i++) {
        _atoms_ratio[i] = ratio[i];
    }
//  fixme  this->_mass = mass;
//  fixme  this->_mass_factor = 1 / sqrt(mass);
    return *this;
}

WorldBuilder &WorldBuilder::setBoxSize(int64_t box_x, int64_t box_y, int64_t box_z) {
    this->box_x = box_x;
    this->box_y = box_y;
    this->box_z = box_z;
    return *this;
}

void WorldBuilder::build() {
    if (_p_atom == nullptr) {
        // todo return error.
    }
    if (_p_domain == nullptr) {
        // todo return error
    }

    createPhaseSpace();

    double p[4] = {0.0, 0.0, 0.0, 0.0}; // index 0-2: mv in 3 d; index 3: mass total at local.
    double _vcm[4]; // same as @var p, but global.
    vcm(p);
    MPI_Allreduce(p, _vcm, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    _type_atom_count n_atoms_global = 2 * (unsigned long) box_x * (unsigned long) box_y * (unsigned long) box_z;
//    double mass_total = n_atoms_global * atom_type::getAtomMass(atom_type::Fe); // todo multiple type.
    double &mass_total = _vcm[3];

    if (mass_total > 0.0) {
        _vcm[0] /= n_atoms_global; // the momentum to be cut off.
        _vcm[1] /= n_atoms_global;
        _vcm[2] /= n_atoms_global;
    }

    zeroMomentum(_vcm);
#ifdef DEV_MODE
    vcm(p);
    kiwi::logs::d("momentum", "momentum:{0} {1} {2}\n", p[0], p[1], p[2]);
#endif

    // double scalar, tfactor;
    // double t_test, scalar_test;
    // t_test = t;
    double scalar = computeScalar(n_atoms_global);

    //t_test *= tfactor;
    //MPI_Allreduce(&t_test,&scalar_test,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    /**
     * \sum m_i(v_i)^2 = 3*nkT  => rescale_factor^2 = T = \sum m_i(v_i)^2 / 3nk
     * thus: T / T_set =  \sum m_i(v_i)^2 / \sum m_i(v'_i)^2
     * then:
     */
    double rescale_factor = sqrt(tset / scalar);
    rescale(rescale_factor); // todo
}

void WorldBuilder::createPhaseSpace() {
    unsigned long id_pre = (unsigned long) box_x * box_y * _p_domain->getSubBoxLatticeCoordLower(2)
                           + (unsigned long) _p_domain->getSubBoxLatticeCoordLower(1) *
                             box_x * _p_domain->getSubBoxLatticeSize(2)
                           + (unsigned long) _p_domain->getSubBoxLatticeCoordLower(0) *
                             _p_domain->getSubBoxLatticeSize(1) * _p_domain->getSubBoxLatticeSize(2);
    /*for(int i = 0; i < id_pre; i++){
        uniform();
        uniform();
        uniform();
    }*/
    _type_atom_mass mass = 0;
    for (int k = 0; k < _p_domain->getSubBoxLatticeSize(2); k++) {
        for (int j = 0; j < _p_domain->getSubBoxLatticeSize(1); j++) {
            for (int i = 0; i < _p_domain->getSubBoxLatticeSize(0); i++) {
                AtomElement &atom_ = _p_atom->getAtomList()->getAtomEleBySubBoxIndex(i, j, k);
//                kk = _p_atom->IndexOf3DIndex(i, j, k); // todo
                atom_.id = ++id_pre;
                atom_._tp = randomAtomsType(); // set random atom type.
                mass = atom_type::getAtomMass(atom_._tp); // get atom mass of this kind of atom.
                // atoms[kk].x[0] = (_p_domain->getSubBoxLatticeCoordLower(0) + (i - xstart)) * (_lattice_const / 2);
                atom_.x[0] = (_p_domain->getSubBoxLatticeCoordLower(0) + i) * 0.5 * (_lattice_const);
                atom_.x[1] = (_p_domain->getSubBoxLatticeCoordLower(1) + j) * _lattice_const +
                             (i % 2) * (_lattice_const / 2);
                atom_.x[2] = (_p_domain->getSubBoxLatticeCoordLower(2) + k) * _lattice_const +
                             (i % 2) * (_lattice_const / 2);
                atom_.v[0] = (uniform() - 0.5) / mass;
                atom_.v[1] = (uniform() - 0.5) / mass;
                atom_.v[2] = (uniform() - 0.5) / mass;
            }
        }
    }
}

/**
 * To achieve the goal of zero momentum,
 * we denote the total momentum of the system as vcm, the total mass of all atoms in system as M.
 * And we assume the total count of atoms is N.
 *
 * We cutoff the momentum of each atom by vcm/N, that is $ v_i' = v_i - vcm/(N * m_i) $.
 * in which, v_i is the velocity before cutting of, v_i' is the velocity after cutting of.
 * Because, \sum_{i=1}^{N} v_i * m_i = \sum_{i=1}^{N} vcm/N = vcm.
 * Thus, \sum_{i=1}^{N} v_i' * m_i  =0.
 */
void WorldBuilder::zeroMomentum(double *vcm) {
//    int xstart = _p_domain->getSubBoxLatticeCoordLower(0) - _p_domain->getGhostLatticeCoordLower(0);
//    int ystart = _p_domain->getSubBoxLatticeCoordLower(1) - _p_domain->getGhostLatticeCoordLower(1);
//    int zstart = _p_domain->getSubBoxLatticeCoordLower(2) - _p_domain->getGhostLatticeCoordLower(2);
//    long kk; // todo unsigned
    _type_atom_mass mass;
    for (int k = 0; k < _p_domain->getSubBoxLatticeSize(2); k++) {
        for (int j = 0; j < _p_domain->getSubBoxLatticeSize(1); j++) {
            for (int i = 0; i < _p_domain->getSubBoxLatticeSize(0); i++) {
                AtomElement &atom_ = _p_atom->getAtomList()->getAtomEleBySubBoxIndex(i, j, k);
//                kk = _p_atom->IndexOf3DIndex(i, j, k);
                mass = atom_type::getAtomMass(atom_._tp);
                atom_.v[0] -= (vcm[0] / mass);
                atom_.v[1] -= (vcm[1] / mass);
                atom_.v[2] -= (vcm[2] / mass);
            }
        }
    }
}

double WorldBuilder::computeScalar(_type_atom_count n_atoms) {
    double t = 0.0;
//    int xstart = _p_domain->getSubBoxLatticeCoordLower(0) - _p_domain->getGhostLatticeCoordLower(0);
//    int ystart = _p_domain->getSubBoxLatticeCoordLower(1) - _p_domain->getGhostLatticeCoordLower(1);
//    int zstart = _p_domain->getSubBoxLatticeCoordLower(2) - _p_domain->getGhostLatticeCoordLower(2);
//    long kk; // todo unsigned
    for (int k = 0; k < _p_domain->getSubBoxLatticeSize(2); k++) {
        for (int j = 0; j < _p_domain->getSubBoxLatticeSize(1); j++) {
            for (int i = 0; i < _p_domain->getSubBoxLatticeSize(0); i++) {
                AtomElement &atom_ = _p_atom->getAtomList()->getAtomEleBySubBoxIndex(i, j, k);
//                kk = _p_atom->IndexOf3DIndex(i, j, k);
                t += (atom_.v[0] * atom_.v[0] +
                      atom_.v[1] * atom_.v[1] +
                      atom_.v[2] * atom_.v[2]) * atom_type::getAtomMass(atom_._tp); // fixme
            }
        }
    }

    double t_global;
    MPI_Allreduce(&t, &t_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    t_global *= dofCompute(n_atoms);
    return t_global;
}

void WorldBuilder::rescale(double rescale_factor) {
//    int xstart = _p_domain->getSubBoxLatticeCoordLower(0) - _p_domain->getGhostLatticeCoordLower(0);
//    int ystart = _p_domain->getSubBoxLatticeCoordLower(1) - _p_domain->getGhostLatticeCoordLower(1);
//    int zstart = _p_domain->getSubBoxLatticeCoordLower(2) - _p_domain->getGhostLatticeCoordLower(2);
//    long kk;
    for (int k = 0; k < _p_domain->getSubBoxLatticeSize(2); k++) {
        for (int j = 0; j < _p_domain->getSubBoxLatticeSize(1); j++) {
            for (int i = 0; i < _p_domain->getSubBoxLatticeSize(0); i++) {
                AtomElement &atom_ = _p_atom->getAtomList()->getAtomEleBySubBoxIndex(i, j, k);
//                kk = _p_atom->IndexOf3DIndex(i, j, k);
                atom_.v[0] *= rescale_factor;
                atom_.v[1] *= rescale_factor;
                atom_.v[2] *= rescale_factor;
            }
        }
    }
}

double WorldBuilder::dofCompute(unsigned long n_atoms) {
    unsigned long dof = 3 * n_atoms;
    dof -= 3; // fixme, why?
    return mvv2e / (dof * BOLTZ);
}

double WorldBuilder::uniform() {
//#ifdef DEV_MODE
//    return 1;
//#else
    int k = _random_seed / IQ;
    _random_seed = IA * (_random_seed - k * IQ) - IR * k;
    if (_random_seed < 0) {
        _random_seed += IM;
    }
    double ans = AM * _random_seed;
    return ans;
//#endif
}

void WorldBuilder::vcm(double p[DIMENSION + 1]) {
    _type_atom_mass mass_one = 0.0;
//    int xstart = _p_domain->getSubBoxLatticeCoordLower(0) - _p_domain->getGhostLatticeCoordLower(0);
//    int ystart = _p_domain->getSubBoxLatticeCoordLower(1) - _p_domain->getGhostLatticeCoordLower(1);
//    int zstart = _p_domain->getSubBoxLatticeCoordLower(2) - _p_domain->getGhostLatticeCoordLower(2);
    // reset p.
    for (int i = 0; i < DIMENSION; i++) {
        p[i] = 0;
    }
    p[DIMENSION] = 0;

//    long kk = 0;
    for (int k = 0; k < _p_domain->getSubBoxLatticeSize(2); k++) {
        for (int j = 0; j < _p_domain->getSubBoxLatticeSize(1); j++) {
            for (int i = 0; i < _p_domain->getSubBoxLatticeSize(0); i++) {
                AtomElement &atom_ = _p_atom->getAtomList()->getAtomEleBySubBoxIndex(i, j, k);
//                kk = _p_atom->IndexOf3DIndex(i, j, k);
                mass_one = atom_type::getAtomMass(atom_._tp);
                p[0] += atom_.v[0] * mass_one;
                p[1] += atom_.v[1] * mass_one;
                p[2] += atom_.v[2] * mass_one;
                p[3] += mass_one; // all mass.
//                kk++;
            }
        }
    }
}

atom_type::atom_type WorldBuilder::randomAtomsType() {
    int ratio_total = 0;
    for (int i = 0; i < atom_type::num_atom_types; i++) {
        ratio_total += _atoms_ratio[i];
    }
#ifdef DEV_MODE
    int rand_ = rand() % ratio_total;
#else
    int rand_ = rand() % ratio_total; // todo srank. Rand() has limited randomness; use C++ lib instead.
#endif
//    return atom_type::getAtomTypeByOrder();
    int rank_local = 0;
    for (int i = 0; i < atom_type::num_atom_types; i++) {
        rank_local += _atoms_ratio[i];
        if (rand_ < rank_local) {
            return atom_type::getAtomTypeByNum(i);
        }
    }
//    return atom_type::Fe;
}
