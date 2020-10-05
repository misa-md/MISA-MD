//
// Created by genshen on 2018-05-05.
// based on create_atom.h/create_atom.cpp created by baihe back to 2016-01-18.
//
//

#ifndef MISA_MD_WORLD_BUILDER_H
#define MISA_MD_WORLD_BUILDER_H

#include <comm/domain/domain.h>
#include "atom.h"
#include "types/atom_types.h"

// todo documents
class WorldBuilder {
public:
    WorldBuilder();

    WorldBuilder &setDomain(comm::BccDomain *p_domain);

    WorldBuilder &setAtomsContainer(AtomSet *p_atom);

    /**
     * initialize random seed
     * (note: in md, random is only used to create atoms)
     * @param seek random seed
     * @return reference of builder.
     */
    WorldBuilder &setRandomSeed(int seek);

    /**
     * set temperature for creation.
     * @param t_set
     * @return
     */
    WorldBuilder &setTset(double t_set);

    WorldBuilder &setLatticeConst(double lattice_const);

    /**
     * set the ratio of alloy(e.g. Fe-Cu-Ni alloy)
     * @param mass
     * @return
     */
    WorldBuilder &setAlloyRatio(const int ratio[atom_type::num_atom_types]);

    WorldBuilder &setBoxSize(int64_t box_x, int64_t box_y, int64_t box_z);

    /**
     * set position and random velocity of all atoms in zero-momentum condition
     */
    void build();

    /**
     * calculate the total momentum of atom in this sub-box, returned by array p.
     * @param p for return value. p[0] to p[DIMENSION -1] is the momentum at each dimension;
     *          p[DIMENSION] is the total mass of atoms in this box.
     */
    void vcm(double p[DIMENSION + 1]);

protected:
    atom_type::atom_type randomAtomsType();

private:
    comm::BccDomain *_p_domain;
    AtomSet *_p_atom;

    int64_t box_x = 0, box_y = 0, box_z = 0; // todo re type
    double tset;
    double _lattice_const;
    int _atoms_ratio[atom_type::num_atom_types];
//    fixme double _mass, _mass_factor;
//    double _mass, _mass_factor;

    void createPhaseSpace();

    /**
     * change velocity of each atom to make the total momentum of the system equals to zero.
     * @param vcm the total momentum of system at each dimension.
     */
    void zeroMomentum(double *vcm);

};


#endif //MISA_MD_WORLD_BUILDER_H
