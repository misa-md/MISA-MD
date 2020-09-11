//
// Created by genshen on 2019-03-17.
//

#ifndef MISA_MD_LAT_PARTICLE_PACKER_H
#define MISA_MD_LAT_PARTICLE_PACKER_H

#include <vector>
#include <comm/packer.h>
#include <comm/domain/domain.h>
#include "atom/atom_list.h"
#include "lat_particle_data.h"


/**
 * lattice particle packer for exchange lattice atoms with neighbours.
 */
class LatParticlePacker : public comm::Packer<LatParticleData> {
public:
    LatParticlePacker(const comm::BccDomain &domain, AtomList &atom_list,
                      std::vector<std::vector<_type_atom_id>> &send_list,
                      std::vector<std::vector<_type_atom_id>> &receive_list);

    const unsigned long sendLength(const int dimension, const int direction) override;

    void setOffset(double offset[DIMENSION],
                   const int dimension, const int direction);

protected:
    const comm::BccDomain &domain;
    AtomList &atom_list;
    std::vector<std::vector<_type_atom_id>> &send_list;
    std::vector<std::vector<_type_atom_id>> &receive_list;
};

class LatPackerFirst : public LatParticlePacker {
public:
    LatPackerFirst(const comm::BccDomain &domain, AtomList &atom_list,
                   std::vector<std::vector<_type_atom_id>> &send_list,
                   std::vector<std::vector<_type_atom_id>> &receive_list)
            : LatParticlePacker(domain, atom_list, send_list, receive_list) {};

    void onSend(LatParticleData buffer[], const unsigned long send_len,
                const int dimension, const int direction) override;

    /**
     * see {@path figure_pack_area.txt} for detail.
     * In the figure, at x dimension, the packed atoms is in box A (front and back);
     * at y dimension, the packed atoms is in box B (left and right);
     * at z dimension, the packed atoms is in box C (up and bottom);
     * the volume: V(A) < V(B) < V(C).
     * @param buffer data buffer received.
     * @param receive_len the count of data to be packed.
     * @param dimension dimension {0,1,2}
     * @param direction direction 0,1; left (0) or right, up(1) or bottom(0), front(1) or back(0)
     */
    void onReceive(LatParticleData buffer[], const unsigned long receive_len,
                   const int dimension, const int direction) override;
};

class LatPacker : public LatParticlePacker {
public:
    LatPacker(const comm::BccDomain &domain, AtomList &atom_list,
              std::vector<std::vector<_type_atom_id>> &send_list,
              std::vector<std::vector<_type_atom_id>> &receive_list)
            : LatParticlePacker(domain, atom_list, send_list, receive_list) {};

    void onSend(LatParticleData buffer[], const unsigned long send_len,
                const int dimension, const int direction) override;

    void onReceive(LatParticleData buffer[], const unsigned long receive_len,
                   const int dimension, const int direction) override;
};


#endif //MISA_MD_LAT_PARTICLE_PACKER_H
