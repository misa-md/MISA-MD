//
// Created by baihe back to 2016-01-06.
//

#ifndef CRYSTAL_MD_ATOM_H
#define CRYSTAL_MD_ATOM_H

class DomainDecomposition;

#include <cstdio>
#include <vector>
#include <io/io_writer.h>
#include "domain.h"
#include "eam.h"
#include "particledata.h"
#include "lat_particle_data.h"


#define COORDINATE_ATOM_OUT_BOX (-100)

using namespace std;

class atom {
public :
    atom(DomainDecomposition *domain, double latticeconst,
         double cutoffRadius, int seed);

    ~atom();

    void addatom(unsigned long id, double rx, double ry, double rz, double vx, double vy, double vz);

    int decide();

    void clearForce();

    void computeEam(eam *pot, DomainDecomposition *_domaindecomposition, double &comm);

    int getinteridsendsize();

    void computefirst(double dtInv2m, double dt);

    void computesecond(double dtInv2m);

    void getatomx(int direction, vector<vector<int>> &sendlist);

    void getatomy(int direction, vector<vector<int>> &sendlist);

    void getatomz(int direction, vector<vector<int>> &sendlist);

    void getIntertosend(int d, int direction, double ghostlengh, vector<int> &sendlist);

    int getintersendnum(int dimension, int direction);

    void pack_intersend(particledata *buf);

    void unpack_interrecv(int d, int n, particledata *buf);

    void pack_bordersend(int dimension, int n, vector<int> &sendlist, LatParticleData *buf, double shift);

    void unpack_borderrecv(int n, LatParticleData *buf, vector<int> &recvlist);

    void pack_send(int dimension, int n, vector<int> &sendlist, LatParticleData *buf, double shift);

    void unpack_recvfirst(int d, int direction, int n, LatParticleData *buf, vector<vector<int> > &recvlist);

    void unpack_recv(int d, int direction, int n, LatParticleData *buf, vector<vector<int>> &recvlist);

    void pack_rho(int n, vector<int> &recvlist, double *buf);

    void unpack_rho(int d, int direction, double *buf, vector<vector<int>> &sendlist);

    void pack_df(vector<int> &sendlist, vector<int> &intersendlist, double *buf);

    void unpack_df(int n, double *buf, vector<int> &recvlist, vector<int> &interrecvlist);

    void pack_force(int n, vector<int> &recvlist, double *buf);

    void unpack_force(int d, int direction, double *buf, vector<vector<int>> &sendlist);

    void printAtoms(int rank, int outMode, kiwi::IOWriter *writer);

    void setv(int lat[4], double collision_v[3]);

    int getnlocalatom();

    void createphasespace(double factor, int box_x, int box_y, int box_z);

    void vcm(double mass, double masstotal, double *p);

    void zero_momentum(double *vcm);

    double compute_scalar(double mass);

    void rescale(double scalar, double t_set);

    void print_force();

private:
    void calculateNeighbourIndices();

    long int IndexOf3DIndex(long int xIndex, long int yIndex, long int zIndex) const;

    double uniform();

    DomainDecomposition *_domain;
    int numberoflattice;

    double _cutoffRadius;
    int _cutlattice;
    double _latticeconst;
    int _seed;

    int nlocalx, nlocaly, nlocalz; // 本地box内晶格数
    int nghostx, nghosty, nghostz; // ghost区域+local区域内晶格数

    int lolocalx, lolocaly, lolocalz;
    int loghostx, loghosty, loghostz; // 本地对应的全局晶格坐标

    vector<long int> NeighbourOffsets; // 邻居粒子偏移量

    //晶格点原子用数组存储其信息
    unsigned long *id;
    int *type;
    double *x, *v, *f, *rho, *df;

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
