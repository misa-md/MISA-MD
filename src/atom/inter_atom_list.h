//
// Created by genshen on 2018-05-19.
//

#ifndef CRYSTAL_MD_INTER_ATOM_LIST_H
#define CRYSTAL_MD_INTER_ATOM_LIST_H

#include <vector>
#include <list>
#include "domain.h"
#include "atom_element.h"
#include "../pack/particledata.h"
#include "../pack/lat_particle_data.h"
#include "../types/pre_define.h"
#include "../types/atom_types.h"

typedef std::list<AtomElement> _type_inter_list;

/**
 * storing inter atoms
 */
class InterAtomList {
    friend class atom;

public:

    InterAtomList();

    void appendInter(_type_atom_id atom_id);

    /**
     * add an inter atom to @var inter_list.
     * The atom data will be copied into the list.
     * @param atom reference of the atom.
     */
    void addInterAtom(AtomElement &atom);

    inline size_t nLocalInter() {
        return nlocalinter;
    }

    void exchangeInter(Domain *p_domain);

    void borderInter(Domain *p_domain);

    _type_inter_list inter_list;
    _type_inter_list inter_ghost_list;
    size_t nlocalinter; // 本地间隙原子数
    size_t nghostinter; // ghost间隙原子数

    /**
     * pointer of element in atom_list (pointer of {@class AtomElement}).
     * // todo use avl tree.
     * // todo use pointer.
     */
private:
    std::vector<std::vector<int> > intersendlist; // todo make it temp variable.
    std::vector<std::vector<int> > interrecvlist; // todo make it temp variable.
    std::vector<unsigned long> interbuf; // todo make it temp variable.

    void pack_intersend(std::vector<unsigned long> interbuf, particledata *buf);

    void unpack_interrecv(int d, int n,
                          const double lower[DIMENSION], // p_domain->getMeasuredSubBoxLowerBounding(d)
                          const double upper[DIMENSION], // p_domain->getMeasuredSubBoxUpperBounding(d)
                          particledata *buf);

    void pack_bordersend(int dimension, int n, std::vector<int> &sendlist, LatParticleData *buf, double shift);

    void unpack_borderrecv(int n, const double lower[DIMENSION], // p_domain->getMeasuredGhostLowerBounding(d)
                           const double upper[DIMENSION], // p_domain->getMeasuredGhostUpperBounding(d)
                           LatParticleData *buf, std::vector<int> &recvlist);

    unsigned long getinteridsendsize();

    /**
     * If some inter atoms get into ghost area of neighbour processors(they are still in local box.),
     * those atoms should be send to neighbour processors
     * (neighbour processors will save those atoms as ghost intel atoms).
     * We call those atoms as "neighbour ghost intel atom".
     *
     * @brief This method will record those atoms.
     * @param p_domain pointer of simulation domain
     * @param d dimension 0,1,2 of 3d. @param d values = {0,1,2}
     * @param direction direction of LOW or HIGH. One direction has 2 direction(such as up and down, back and front, left and right).
     * @param ghostlengh the measured length of ghost area.
     * @param sendlist the atoms to be send will be saved in this data.
     */
    void getIntertosend(Domain *p_domain, int d, int direction, double ghostlengh, std::vector<int> &sendlist);

    /**
     * If some inter atoms get out of box, those atom is no more in current dox of current processor.
     * they should be send to corresponding neighbour processors.
     *
     * @param p_domain pointer of simulation domain
     * @param dimension dimension of 3d. @param d values = {0,1,2}
     * @param direction direction of LOW or HIGH.
     * @return the number of out-of-box inter atoms.
     */
    unsigned long getintersendnum(Domain *p_domain, int dimension, int direction);
};


#endif //CRYSTAL_MD_INTER_ATOM_LIST_H
