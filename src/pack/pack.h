//
// Created by genshen on 5/17/18.
//

#ifndef CRYSTAL_MD_PACK_H
#define CRYSTAL_MD_PACK_H

#include <vector>
#include "../config/pre_define.h"
#include "../atom/atom_list.h"
#include "lat_particle_data.h"
#include "particledata.h"
#include "../atom/inter_atom_list.h"

namespace pack {

    void pack_intersend(InterAtomList *inter,
                        std::vector<unsigned long> interbuf, particledata *buf);

    void unpack_interrecv(int d, int n, InterAtomList *inter,
                          double lower[DIMENSION], // p_domain->getMeasuredSubBoxLowerBounding(d)
                          double upper[DIMENSION], // p_domain->getMeasuredSubBoxUpperBounding(d)
                          particledata *buf);

    void pack_bordersend(int dimension, int n, InterAtomList *inter,
                         std::vector<int> &sendlist, LatParticleData *buf, double shift);

    void unpack_borderrecv(int n, InterAtomList *inter,
                           double lower[DIMENSION], // p_domain->getMeasuredGhostLowerBounding(d)
                           double upper[DIMENSION], // p_domain->getMeasuredGhostUpperBounding(d)
                           LatParticleData *buf, std::vector<int> &recvlist);

    /**
     * package ghost atom to its neighbors processors
     * @param dimension 0,1,2. which refers to x,y,z dimension.
     * @param n the atoms count to be packed.
     * @param sendlist id list of atoms to be packed.
     * @param buf buffer to store packed ghost atoms data (e.g. atom type and atom location).
     * @param shift coordinate offset used for periodic boundary.
     * e.g: it will add [global box length] to the coordinate of ghost atoms at leftmost sub-box to fit periodic boundary.
     */
    void pack_send(const int dimension, const int n, const double shift[DIMENSION], AtomList &atom_list,
                   LatParticleData *buf, std::vector<_type_atom_id> &sendlist);

    /**
     * see {@path figure_pack_area.txt} for detail.
     * In the figure, at x dimension, the packed atoms is in box A (front and back);
     * at y dimension, the packed atoms is in box B (left and right);
     * at z dimension, the packed atoms is in box C (up and bottom);
     * the volume: V(A) < V(B) < V(C).
     * @param d dimension (0,1,2)
     * @param direction direction 0,1; left (0) or right, up(1) or bottom(0), front(1) or back(0)
     * @param n the count of data to be packed.
     * @param atom_list
     * @param buf
     * @param recvlist
     */
    void unpack_recvfirst(int d, int direction, int n, AtomList &atom_list,
                          _type_lattice_size ghost[DIMENSION], //p_domain->getGhostLatticeSize(d)
                          _type_lattice_size box[DIMENSION], //p_domain->getSubBoxLatticeSize(d)
                          _type_lattice_size ext[DIMENSION], //p_domain->getGhostExtLatticeSize(d)
                          LatParticleData *buf, std::vector<std::vector<_type_atom_id> > &recvlist);

    void unpack_recv(int d, int direction, int n, AtomList &atom_list,
                     LatParticleData *buf, std::vector<std::vector<_type_atom_id>> &recvlist);

    void pack_rho(int n, AtomList &atom_list, double *buf, std::vector<_type_atom_id> &recvlist);

    void unpack_rho(int d, int direction, AtomList &atom_list,
                    double *buf, std::vector<std::vector<_type_atom_id>> &sendlist);

    void pack_df(AtomList &atom_list, double *buf, InterAtomList *inter,
                 std::vector<_type_atom_id> &sendlist, std::vector<int> &intersendlist);

    void unpack_df(int n, AtomList &atom_list, double *buf, InterAtomList *inter,
                   std::vector<_type_atom_id> &recvlist, std::vector<int> &interrecvlist);

    void pack_force(int n, AtomList &atom_list, double *buf, std::vector<_type_atom_id> &recvlist);

    void unpack_force(int d, int direction, AtomList &atom_list, double *buf,
                      std::vector<std::vector<_type_atom_id> > &sendlist);
};


#endif //CRYSTAL_MD_PACK_H
