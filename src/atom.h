//
// Created by baihe back to 2016-01-06.
//

#ifndef CRYSTAL_MD_ATOM_H
#define CRYSTAL_MD_ATOM_H

#include <cstdio>
#include <vector>
#include <eam.h>

#include "domain/domain.h"
#include "atom/atom_element.h"
#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"
#include "pack/particledata.h"
#include "pack/lat_particle_data.h"

class Domain; // todo remove.

class atom {
public :
    friend class Domain;

    atom(Domain *domain, double latticeconst,
         double cutoffRadiusFactor, int seed);

    ~atom();

    /**
     * used in read creating mode.
     */
    void addAtom(unsigned long id, double rx, double ry, double rz, double vx, double vy, double vz);

    /**
     * compute the index offset of neighbour atoms.
     */
    void calculateNeighbourIndices();

    int decide();

    void clearForce();

    void computeEam(eam *pot, Domain *domain, double &comm);

    /**
     * set velocity of a atom whose position is specified by array @param lat
     * This atom is called PKA (Primary Knock-on Atom).
     * We first compute the the velocity in each direction from the energy, and the set velocity of the atom.
     * The unit of energy is ev, the unit of velocity is 100m/s.
     *
     * @param lat the position of pka.
     * @param direction the vector of velocity to be set.
     * @param energy the energy of PKA.
     */
    void setv(int lat[4], double direction[3], double energy);

    int getnlocalatom();

    void print_force();

    inline AtomList *getAtomList() {
        return atom_list;
    }

    inline AtomList &getAtomListRef() {
        return *atom_list;
    }

    inline InterAtomList *getInterList() {
        return inter_atom_list;
    }

    inline int getCutLattice() {
        return _cutlattice;
    }

    void sendForce();

private:

    Domain *p_domain;
    int numberoflattice;

    double _cutoffRadius;
    int _cutlattice;
    double _latticeconst;
    int _seed;

    std::vector<long int> NeighbourOffsets; // 邻居粒子偏移量 // todo use offset in x,y,z dimension

    AtomList *atom_list;
    InterAtomList *inter_atom_list;

    void sendrho();

    void sendDfEmbed();

};

#endif // CRYSTAL_MD_ATOM_H
