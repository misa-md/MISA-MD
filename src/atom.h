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
#include "particledata.h"
#include "lat_particle_data.h"
#include "atom_list.h"

class Domain; // todo remove.

class atom {
public :

    friend class AtomDump;

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

    void pack_intersend(particledata *buf);

    void unpack_interrecv(int d, int n, particledata *buf);

    void pack_bordersend(int dimension, int n, vector<int> &sendlist, LatParticleData *buf, double shift);

    void unpack_borderrecv(int n, LatParticleData *buf, vector<int> &recvlist);

    /**
     * package ghost atom to its neighbors processors
     * @param dimension 0,1,2. which refers to x,y,z dimension.
     * @param n the atoms count to be packed.
     * @param sendlist id list of atoms to be packed.
     * @param buf buffer to store packed ghost atoms data (e.g. atom type and atom location).
     * @param shift coordinate offset used for periodic boundary.
     * e.g: it will add [global box length] to the coordinate of ghost atoms at leftmost sub-box to fit periodic boundary.
     */
    void pack_send(int dimension, int n, vector<_type_atom_id> &sendlist, LatParticleData *buf, double shift[DIMENSION]);

    void unpack_recvfirst(int d, int direction, int n, LatParticleData *buf, vector<vector<_type_atom_id> > &recvlist);

    void unpack_recv(int d, int direction, int n, LatParticleData *buf, vector<vector<_type_atom_id>> &recvlist);

    void pack_rho(int n, vector<_type_atom_id> &recvlist, double *buf);

    void unpack_rho(int d, int direction, double *buf, vector<vector<_type_atom_id>> &sendlist);

    void pack_df(vector<_type_atom_id> &sendlist, vector<int> &intersendlist, double *buf);

    void unpack_df(int n, double *buf, vector<_type_atom_id> &recvlist, vector<int> &interrecvlist);

    void pack_force(int n, vector<_type_atom_id> &recvlist, double *buf);

    void unpack_force(int d, int direction, double *buf, vector<vector<_type_atom_id>> &sendlist);

    void setv(int lat[4], double collision_v[3]);

    int getnlocalatom();

    void print_force();

    AtomList *getAtomList() {
        return atom_list;
    }

private:

    long IndexOf3DIndex(long int xIndex, long int yIndex, long int zIndex) const;

//    double uniform();

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
    vector<unsigned long> idinter;
    vector<int> typeinter;
    vector<vector<double>> xinter; // 间隙原子坐标
    vector<vector<double>> vinter; // 间隙原子速度
    vector<vector<double>> finter; // 间隙原子力
    vector<double> rhointer;
    vector<double> dfinter;
    int nlocalinter, nghostinter; // 本地间隙原子数和ghost间隙原子数

    vector<unsigned long> interbuf;
};

#endif // CRYSTAL_MD_ATOM_H
