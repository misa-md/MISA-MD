//
// Created by genshen on 2018-05-05.
// based on create_atom.h/create_atom.cpp created by baihe back to 2016-01-18.
//
//

#ifndef CRYSTALMD_WORLD_BUILDER_H
#define CRYSTALMD_WORLD_BUILDER_H

#include "atom.h"
#include "atom_types.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

// todo documents
class WorldBuilder {
public:
    WorldBuilder();

    WorldBuilder &setDomain(Domain *p_domain);

    WorldBuilder &setAtomsContainer(atom *p_atom);

    WorldBuilder &setRandomSeed(int seek);

    WorldBuilder &setTset(double tset);

    WorldBuilder &setLatticeConst(double lattice_const);

    /**
     * set the ratio of alloy(e.g. Fe-Cu-Ni alloy)
     * @param mass
     * @return
     */
    WorldBuilder &setAlloyRatio(int ratio[atom_type::num_atom_types]);

    WorldBuilder &setBoxSize(int box_x, int box_y, int box_z);

    void build();

private:
    Domain *_p_domain;
    atom *_p_atom;

    int _random_seed; // random seed for creating atoms.
    int box_x = 0, box_y = 0, box_z = 0;
    double tset;
    double _lattice_const;
    int _atoms_ratio[atom_type::num_atom_types];
//    fixme double _mass, _mass_factor;
    double _mass, _mass_factor;

    double dofCompute(unsigned long natom);

    void createPhaseSpace();

    void vcm(double mass, double masstotal, double *p);

    double uniform();

    /**
     * change velocity of each atom to make the total momentum of the system equals to zero.
     * @param vcm the total momentum of system at each dimension.
     */
    void zeroMomentum(double *vcm);

    double computeScalar();

    void rescale(double scalar);
};


#endif //CRYSTALMD_WORLD_BUILDER_H
