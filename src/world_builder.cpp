//
// Created by genshen on 5/5/18.
//

#include <cmath>
#include "world_builder.h"

WorldBuilder::WorldBuilder() : _random_seed(1024), tset(0), _mass(1), _mass_factor(1),
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

WorldBuilder &WorldBuilder::setMass(double mass) {
    this->_mass = mass;
    this->_mass_factor = 1 / sqrt(mass);
    return *this;
}

WorldBuilder &WorldBuilder::setBoxSize(int box_x, int box_y, int box_z) {
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

    unsigned long n_atoms = 2 * (unsigned long) box_x * (unsigned long) box_y * (unsigned long) box_z;
    double mass_total = n_atoms * _mass; // todo multiple type.

    double p[3] = {0.0, 0.0, 0.0}, _vcm[3];
    vcm(_mass, mass_total, p);
    MPI_Allreduce(p, _vcm, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (mass_total > 0.0) {
        _vcm[0] /= mass_total;
        _vcm[1] /= mass_total;
        _vcm[2] /= mass_total;
    }

    zeroMomentum(_vcm);

    double scalar, tfactor;
    //double t_test, scalar_test;
    //t_test = t;
    double t = computeScalar();
    MPI_Allreduce(&t, &scalar, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    tfactor = dofCompute(n_atoms);
    scalar *= tfactor;
    //t_test *= tfactor;
    //MPI_Allreduce(&t_test,&scalar_test,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    rescale(scalar); // todo
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
    int xstart = _p_domain->getSubBoxLatticeCoordLower(0) - _p_domain->getGhostLatticeCoordLower(0);
    int ystart = _p_domain->getSubBoxLatticeCoordLower(1) - _p_domain->getGhostLatticeCoordLower(1);
    int zstart = _p_domain->getSubBoxLatticeCoordLower(2) - _p_domain->getGhostLatticeCoordLower(2);
    unsigned long kk;
    for (int k = zstart; k < _p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
        for (int j = ystart; j < _p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
            for (int i = xstart; i < _p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                kk = _p_atom->IndexOf3DIndex(i, j, k); // todo
                _p_atom->atoms[kk].id = ++id_pre;
                // atoms[kk].x[0] = (_p_domain->getSubBoxLatticeCoordLower(0) + (i - xstart)) * (_lattice_const / 2);
                _p_atom->atoms[kk].x[0] = (_p_domain->getGhostLatticeCoordLower(0) + i) * 0.5 * (_lattice_const);
                _p_atom->atoms[kk].x[1] = (_p_domain->getSubBoxLatticeCoordLower(1) + (j - ystart)) * _lattice_const +
                                          (i % 2) * (_lattice_const / 2);
                _p_atom->atoms[kk].x[2] = (_p_domain->getSubBoxLatticeCoordLower(2) + (k - zstart)) * _lattice_const +
                                          (i % 2) * (_lattice_const / 2);
                _p_atom->atoms[kk].v[0] = (uniform() - 0.5) * _mass_factor;
                _p_atom->atoms[kk].v[1] = (uniform() - 0.5) * _mass_factor;
                _p_atom->atoms[kk].v[2] = (uniform() - 0.5) * _mass_factor;
            }
        }
    }
}

/**
 * To achieve the goal of zero momentum,
 * we denote the total momentum of the system as vcm, the total mass of all atoms in system as M.
 * And we assume there only one type of atoms.
 *
 * We cutoff the velocity of each atom by vcm/M, that is $ v_i = v_i - vcm/M $.
 * Because, \sum_{i=1}^{n} v_i * m = \sum_{i=1}^{n} m*vcm/M , in which, vcm = \sum_{i=1}^{n} v_i * m.
 */
void WorldBuilder::zeroMomentum(double *vcm) {
    int xstart = _p_domain->getSubBoxLatticeCoordLower(0) - _p_domain->getGhostLatticeCoordLower(0);
    int ystart = _p_domain->getSubBoxLatticeCoordLower(1) - _p_domain->getGhostLatticeCoordLower(1);
    int zstart = _p_domain->getSubBoxLatticeCoordLower(2) - _p_domain->getGhostLatticeCoordLower(2);
    long kk; // todo unsigned
    for (int k = zstart; k < _p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
        for (int j = ystart; j < _p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
            for (int i = xstart; i < _p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                kk = _p_atom->IndexOf3DIndex(i, j, k);
                _p_atom->atoms[kk].v[0] -= vcm[0];
                _p_atom->atoms[kk].v[1] -= vcm[1];
                _p_atom->atoms[kk].v[2] -= vcm[2];
            }
        }
    }
}

double WorldBuilder::computeScalar() {
    double t = 0.0;
    int xstart = _p_domain->getSubBoxLatticeCoordLower(0) - _p_domain->getGhostLatticeCoordLower(0);
    int ystart = _p_domain->getSubBoxLatticeCoordLower(1) - _p_domain->getGhostLatticeCoordLower(1);
    int zstart = _p_domain->getSubBoxLatticeCoordLower(2) - _p_domain->getGhostLatticeCoordLower(2);
    long kk; // todo unsigned
    for (int k = zstart; k < _p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
        for (int j = ystart; j < _p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
            for (int i = xstart; i < _p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                kk = _p_atom->IndexOf3DIndex(i, j, k);
                t += (_p_atom->atoms[kk].v[0] * _p_atom->atoms[kk].v[0] +
                      _p_atom->atoms[kk].v[1] * _p_atom->atoms[kk].v[1] +
                      _p_atom->atoms[kk].v[2] * _p_atom->atoms[kk].v[2]) * _mass;
            }
        }
    }
    return t;
}

void WorldBuilder::rescale(double scalar) {
    double factor = sqrt(tset / scalar);
    int xstart = _p_domain->getSubBoxLatticeCoordLower(0) - _p_domain->getGhostLatticeCoordLower(0);
    int ystart = _p_domain->getSubBoxLatticeCoordLower(1) - _p_domain->getGhostLatticeCoordLower(1);
    int zstart = _p_domain->getSubBoxLatticeCoordLower(2) - _p_domain->getGhostLatticeCoordLower(2);
    long kk;
    for (int k = zstart; k < _p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
        for (int j = ystart; j < _p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
            for (int i = xstart; i < _p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                kk = _p_atom->IndexOf3DIndex(i, j, k);
                _p_atom->atoms[kk].v[0] *= factor;
                _p_atom->atoms[kk].v[1] *= factor;
                _p_atom->atoms[kk].v[2] *= factor;
            }
        }
    }
}

double WorldBuilder::dofCompute(unsigned long natom) {
    unsigned long dof = 3 * natom;
    dof -= 3;
    return mvv2e / (dof * BOLTZ);
}

double WorldBuilder::uniform() {
    int k = _random_seed / IQ;
    _random_seed = IA * (_random_seed - k * IQ) - IR * k;
    if (_random_seed < 0) {
        _random_seed += IM;
    }
    double ans = AM * _random_seed;
    return ans;
}

void WorldBuilder::vcm(double mass, double masstotal, double *p) {
    double massone;
    int xstart = _p_domain->getSubBoxLatticeCoordLower(0) - _p_domain->getGhostLatticeCoordLower(0);
    int ystart = _p_domain->getSubBoxLatticeCoordLower(1) - _p_domain->getGhostLatticeCoordLower(1);
    int zstart = _p_domain->getSubBoxLatticeCoordLower(2) - _p_domain->getGhostLatticeCoordLower(2);
    long kk; // todo unsigned
    for (int k = zstart; k < _p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
        for (int j = ystart; j < _p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
            for (int i = xstart; i < _p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                kk = _p_atom->IndexOf3DIndex(i, j, k);
                massone = mass;
                p[0] += _p_atom->atoms[kk].v[0] * massone;
                p[1] += _p_atom->atoms[kk].v[1] * massone;
                p[2] += _p_atom->atoms[kk].v[2] * massone;
            }
        }
    }
}
