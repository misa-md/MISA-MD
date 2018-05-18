//
// Created by baihe back to 2016-01-06.
//

#ifndef CRYSTAL_MD_ATOM_H
#define CRYSTAL_MD_ATOM_H

#include <cstdio>
#include <vector>
#include <io/io_writer.h>

#include "atom_element.h"
#include "domain.h"
#include "eam.h"
#include "pack/particledata.h"
#include "pack/lat_particle_data.h"
#include "atom_list.h"

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

    void getatomx(int direction, vector<vector<_type_atom_id>> &sendlist);

    void getatomy(int direction, vector<vector<_type_atom_id>> &sendlist);

    void getatomz(int direction, vector<vector<_type_atom_id>> &sendlist);

    void getIntertosend(int d, int direction, double ghostlengh, vector<int> &sendlist);

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

    vector<long int> NeighbourOffsets; // 邻居粒子偏移量

    AtomList *atom_list;

//    unsigned long *id; //
//    int *type;
//    double *x, *v, *f, *rho, *df;
    vector<_type_atom_id > idinter;
    vector<_type_atom_type > typeinter;
    vector<vector<double>> xinter; // 间隙原子坐标
    vector<vector<double>> vinter; // 间隙原子速度
    vector<vector<double>> finter; // 间隙原子力
    vector<double> rhointer;
    vector<double> dfinter;
    int nlocalinter, nghostinter; // 本地间隙原子数和ghost间隙原子数

    vector<unsigned long> interbuf;

    void pack_intersend(particledata *buf);
};

#endif // CRYSTAL_MD_ATOM_H
