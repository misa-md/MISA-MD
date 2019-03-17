//
// Created by genshen on 2019-03-17.
//

#ifndef CRYSTALMD_LAT_PARTICLE_PACKER_H
#define CRYSTALMD_LAT_PARTICLE_PACKER_H

#include <vector>
#include <packer.h>
#include <domain/domain.h>
#include "lat_particle_data.h"

class AtomList;

/**
 * lattice particle packer for exchange lattice atoms with neighbours.
 */
class LatParticlePacker : public Packer<LatParticleData> {
public:
    LatParticlePacker(const comm::Domain &domain, AtomList &atom_list,
                      std::vector<std::vector<_type_atom_id>> &send_list,
                      std::vector<std::vector<_type_atom_id>> &receive_list);

    const unsigned long sendLength(const int dimension, const int direction) override;

    void setOffset(double offset[DIMENSION],
                   const int dimension, const int direction);

protected:
    const comm::Domain &domain;
    AtomList &atom_list;
    std::vector<std::vector<_type_atom_id>> &send_list;
    std::vector<std::vector<_type_atom_id>> &receive_list;
};

class LatPackerFirst : public LatParticlePacker {
public:
    LatPackerFirst(const comm::Domain &domain, AtomList &atom_list,
                   std::vector<std::vector<_type_atom_id>> &send_list,
                   std::vector<std::vector<_type_atom_id>> &receive_list)
            : LatParticlePacker(domain, atom_list, send_list, receive_list) {};

    void onSend(LatParticleData buffer[], const unsigned long send_len,
                const int dimension, const int direction) override;

    void onReceive(LatParticleData buffer[], const unsigned long receive_len,
                   const int dimension, const int direction) override;
};

class LatPacker : public LatParticlePacker {
public:
    LatPacker(const comm::Domain &domain, AtomList &atom_list,
              std::vector<std::vector<_type_atom_id>> &send_list,
              std::vector<std::vector<_type_atom_id>> &receive_list)
            : LatParticlePacker(domain, atom_list, send_list, receive_list) {};

    void onSend(LatParticleData buffer[], const unsigned long send_len,
                const int dimension, const int direction) override;

    void onReceive(LatParticleData buffer[], const unsigned long receive_len,
                   const int dimension, const int direction) override;
};


#endif //CRYSTALMD_LAT_PARTICLE_PACKER_H
