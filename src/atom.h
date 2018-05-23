//
// Created by baihe back to 2016-01-06.
//

#ifndef CRYSTAL_MD_ATOM_H
#define CRYSTAL_MD_ATOM_H

#include <cstdio>
#include <vector>
#include <io/io_writer.h>

#include "domain.h"
#include "atom/atom_element.h"
#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"
#include "pack/particledata.h"
#include "pack/lat_particle_data.h"
#include "potential/eam.h"

class Domain; // todo remove.

class atom {
public :
    friend class AtomDump;
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

    unsigned long getinteridsendsize();

    void computefirst(double dtInv2m, double dt);

    void computesecond(double dtInv2m);

    void getatomx(int direction, std::vector<std::vector<_type_atom_id>> &sendlist);

    void getatomy(int direction, std::vector<std::vector<_type_atom_id>> &sendlist);

    void getatomz(int direction, std::vector<std::vector<_type_atom_id>> &sendlist);

    void getIntertosend(int d, int direction, double ghostlengh, std::vector<int> &sendlist);

    int getintersendnum(int dimension, int direction);

    void setv(int lat[4], double collision_v[3]);

    int getnlocalatom();

    void print_force();

    AtomList *getAtomList() {
        return atom_list;
    }

    AtomList &getAtomListRef() {
        return *atom_list;
    }

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

    std::vector<unsigned long> interbuf;
};

#endif // CRYSTAL_MD_ATOM_H
