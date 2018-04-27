//
// Created by baihe back to 2016-12-22.
//

#ifndef CRYSTAL_MD_DOMAIN_DECOMPOSITION_H
#define CRYSTAL_MD_DOMAIN_DECOMPOSITION_H

#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR

#include <mpi.h>
#include <vector>

class atom;

#include "atom.h"
#include "domain.h"

using namespace std;

#define DIM 3

#define LOWER  0
#define HIGHER 1

/**
 * If N can be decomposed as N = N_x * N_y * N_z, where N, N_x, N_y, N_z are all integer bigger than or equal to 1,
 * then the whole simulation box will be divided into N sub-box with N_z levels in z axis,
 * and in each level, it has  N_x * N_y sub-box.
 * Then, we can bind each processor to a cartesian coordinate (x,y,z) due to the boxes partition,
 * where 0 <= x < N_x, 0 <= z < N_z, 0 <= z < N_z.
 * Last, based on the cartesian coordinate (x,y,z),
 * each processor can get the cartesian coordinate of its contiguous sub-boxes.
 */
class domaindecomposition {
public:
    /**
     * In this construction method, each processor will be bound to a cartesian coordinate.
     *
     * It first divide the simulation box into N pieces(sub-box) (N is the count of all processors).
     * And each processor will be bound to a sub-box, and tagged with a cartesian coordinate(x,y,z).
     */
    domaindecomposition();

    ~domaindecomposition();

    void exchangeInter(atom *_atom, domain *domain);

    void borderInter(atom *_atom, domain *domain);

    void exchangeAtomfirst(atom *_atom, domain *domain);

    void exchangeAtom(atom *_atom, domain *domain);

    void sendrho(atom *_atom);

    void sendDfEmbed(atom *_atom);

    void sendforce(atom *_atom);

    double getBoundingBoxMin(int dimension, domain *domain);

    double getBoundingBoxMax(int dimension, domain *domain);

private:

    MPI_Datatype _mpi_Particle_data;
    MPI_Datatype _mpi_latParticle_data;
    // 每个维度的进程数
    int _gridSize[DIM];
    /**
     * The cartesian coordinate of the sub-box bound to this processor.
     */
    int _coords[DIM];

    // the rank ids of contiguous processors in space.
    int _rank_id_neighbours[DIM][2];

    vector<vector<int> > sendlist;
    vector<vector<int> > recvlist;

    vector<vector<int> > intersendlist;
    vector<vector<int> > interrecvlist;
};

#endif // CRYSTAL_MD_DOMAIN_DECOMPOSITION_H
