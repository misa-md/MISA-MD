//
// Created by genshen on 2018-05-19.
//

#ifndef CRYSTAL_MD_INTER_ATOM_LIST_H
#define CRYSTAL_MD_INTER_ATOM_LIST_H

#include <vector>
#include <list>
#include "domain/domain.h"
#include "atom_element.h"
#include "../pack/particledata.h"
#include "../pack/lat_particle_data.h"
#include "../types/pre_define.h"
#include "../types/atom_types.h"
#include "box.h"

typedef std::list<AtomElement> _type_inter_list;
typedef std::vector<std::vector<AtomElement *> > _type_inter_buf;

/**
 * storing inter atoms
 */
class InterAtomList {
public:
    _type_inter_list inter_list;
    _type_inter_list inter_ghost_list;
    size_t nlocalinter; // 本地间隙原子数
    size_t nghostinter; // ghost间隙原子数

    _type_inter_buf intersendlist;
    _type_inter_buf interrecvlist;

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

    /**
     * pointer of element in atom_list (pointer of {@class AtomElement}).
     * // todo use avl tree.
     * // todo use pointer.
     */
private:
    void pack_bordersend(int dimension, int n, std::vector<AtomElement *> &sendlist,
                         LatParticleData *buf, double shift);

    void unpack_borderrecv(int n, const double lower[DIMENSION], // p_domain->getMeasuredGhostLowerBounding(d)
                           const double upper[DIMENSION], // p_domain->getMeasuredGhostUpperBounding(d)
                           LatParticleData *buf, std::vector<AtomElement *> &recvlist);

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
    void getIntertosend(Domain *p_domain, int d, int direction,
                        double ghostlengh, std::vector<AtomElement *> &sendlist);

    /**
     * count the size of out-of-box inter atoms.
     * @param p_domain  pointer of simulation domain
     * @param n_to_send  the number of out-of-box inter atoms will be stored in each dimension and each direction.
     */
    void countExSendNum(Domain *p_domain, int n_to_send[DIMENSION][2]);

    /**
     * If some inter atoms get out of box, those atom is no more in current dox of current processor.
     * they should be send to corresponding neighbour processors.
     * The out-of-box atoms will be removed from @memberof inter_list in this function.
     *
     * @param p_domain pointer of simulation domain.
     * @param buf the buffer to save the out-of-box atoms.
     * @param dimension dimension of 3d. @param d values = {0,1,2}
     * @param direction direction of LOW or HIGH.
     */
    void packExInterToSend(Domain *p_domain, particledata *buf, int dimension, int direction,
                           box::_type_flag_32 excepted_flag[DIMENSION][2], double offset[DIMENSION]);

    /**
     * unpack exchanged inter atoms data, and save the inter atoms to local inter atom list.
     *
     * @note that we use a delay buffer to recoder the inter atom that is not in this box.
     * those atoms will be written into inter atoms list after exchange finished (delay write).
     * Assume that we write those atoms directly, it will cause some faults:
     * If we save those atoms into inter atom list directly,
     * later communication of other direction or dimension may count those atoms and send them to other processors,
     * but the communication data size (count of atoms to be send to other processors) is precomputed,
     * which can make the data exchange messy.
     *
     * @param buf data buffer
     * @param delay_buffer reference of delay buffer
     * @param n the count of inter atoms in @param buf
     */
    void unpackExInterRecv(Domain *p_domain, particledata *buf,
                           std::list<AtomElement> &delay_buffer,
                           int n);

};


#endif //CRYSTAL_MD_INTER_ATOM_LIST_H
