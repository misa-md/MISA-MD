//
// Created by genshen on 2019-03-17.
//

#ifndef MISA_MD_INTER_BORDER_PACKER_H
#define MISA_MD_INTER_BORDER_PACKER_H


#include <comm/packer.h>
#include <comm/domain/domain.h>

#include "atom/inter_atom_list.h"
#include "particledata.h"
#include "lat_particle_data.h"

/**
 * If some inter atoms get into ghost area of neighbour processors(they are still in local box.),
 * those atoms should be send to neighbour processors
 * (neighbour processors will save those atoms as ghost intel atoms).
 * We call those atoms as "neighbour ghost intel atom".
 *
 */
class InterBorderPacker : public comm::Packer<LatParticleData> {
public:
    explicit InterBorderPacker(const comm::BccDomain &domain, InterAtomList &inter_atom_list);

    /**
     * the atoms to be send will be saved in sendlist..
     * @param dimension dimension 0,1,2 of 3d. @param d values = {0,1,2}
     * @param direction direction of LOW or HIGH. One direction has 2 direction(such as up and down, back and front, left and right).
     * @return the size to be sent.
     */
    const unsigned long sendLength(const int dimension, const int direction) override;

    void onSend(LatParticleData buffer[], const unsigned long send_len,
                const int dimension, const int direction) override;

    void onReceive(LatParticleData buffer[], const unsigned long receive_len,
                   const int dimension, const int direction) override;

private:
    const comm::BccDomain &domain;
    InterAtomList &inter_atom_list;
};


#endif //MISA_MD_INTER_BORDER_PACKER_H
