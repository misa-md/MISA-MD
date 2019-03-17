//
// Created by genshen on 2019-03-17.
//

#ifndef CRYSTALMD_INTER_PARTICLE_PACKER_H
#define CRYSTALMD_INTER_PARTICLE_PACKER_H

#include <list>
#include <packer.h>
#include <domain/domain.h>

#include "atom/atom_element.h"
#include "particledata.h"

class InterAtomList;

class InterParticlePacker : public Packer<particledata> {
public:
    explicit InterParticlePacker(const comm::Domain &domain, InterAtomList &inter_atom_list);

    const unsigned long sendLength(const int dimension, const int direction) override;

    /**
     * If some inter atoms get out of box, those atom is no more in current dox of current processor.
     * they should be send to corresponding neighbour processors.
     * The out-of-box atoms will be removed from @memberof inter_atom_list.inter_list in this function.
     *
     * @param buffer the buffer to save the out-of-box atoms.
     * @param send_len the length of data to be send.
     * @param dimension dimension of 3d. @param d values = {0,1,2}
     * @param direction  direction of LOW or HIGH.
     */
    void onSend(particledata buffer[], const unsigned long send_len,
                const int dimension, const int direction) override;

    /**
     *
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
     * @param n the count of inter atoms in @param buf
     *
     * @param buffer the buffer of received data.
     * @param receive_len the data length received, the count of inter atoms in @param buffer
     * @param dimension dimension of 3d. @param d values = {0,1,2}
     * @param direction  direction of LOW or HIGH.
     */
    void onReceive(particledata buffer[], const unsigned long receive_len,
                   const int dimension, const int direction) override;

    /**
     * do delay write.
     * atoms in @memberof delayed_buffer will be written into inter atoms list in this function.
     */
    void onFinish() override;

protected:
    const comm::Domain &domain;
    InterAtomList &inter_atom_list;
    std::list<AtomElement> delayed_buffer;
    // the size to be send to neighbour processes in each dimension and each direction.
    unsigned long n_to_send[DIMENSION][2] = {{0, 0},
                                             {0, 0},
                                             {0, 0}};
};


#endif //CRYSTALMD_INTER_PARTICLE_PACKER_H
