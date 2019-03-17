//
// Created by genshen on 2019-03-17.
//

#include "lat_particle_packer.h"
#include "pack.h"

LatParticlePacker::LatParticlePacker(const comm::Domain &domain, AtomList &atom_list,
                                     std::vector<std::vector<_type_atom_id>> &send_list,
                                     std::vector<std::vector<_type_atom_id>> &receive_list)
        : domain(domain), atom_list(atom_list),
          send_list(send_list), receive_list(receive_list) {}

const unsigned long LatParticlePacker::sendLength(const int dimension, const int direction) {
    const int index = 2 * dimension + direction;
    return send_list[index].size();
}

void LatParticlePacker::setOffset(double offset[DIMENSION],
                                  const int dimension, const int direction) {
    // 进程在左侧边界
    if (domain.grid_coord_sub_box[dimension] == 0 && direction == LOWER) {
        offset[dimension] = domain.meas_global_length[dimension];
    }
    // 进程在右侧边界
    if (domain.grid_coord_sub_box[dimension] == domain.grid_size[dimension] - 1 && direction == HIGHER) {
        offset[dimension] = -((domain.meas_global_length[dimension]));
    }
}


void LatPackerFirst::onSend(LatParticleData *buffer, const unsigned long send_len,
                            const int dimension, const int direction) {
    // only one periodic boundary appears at one dimension, so the length of array can be 3, not 6.
    double offset[DIMENSION] = {0.0, 0.0, 0.0};
    // 当原子要跨越周期性边界, 原子坐标必须要做出调整
    setOffset(offset, dimension, direction);
    const int index = 2 * dimension + direction;
    pack::pack_send(dimension, send_len, offset, atom_list,
                    buffer, send_list[index]);
}

void LatPackerFirst::onReceive(LatParticleData *buffer, const unsigned long receive_len,
                               const int dimension, const int direction) {
    //将收到的粒子位置信息加到对应存储位置上
    pack::unpack_recvfirst(dimension, direction, receive_len, atom_list,
                           domain.dbx_lattice_size_ghost,
                           domain.dbx_lattice_size_sub_box,
                           domain.dbx_lattice_size_ghost_extended,
                           buffer, receive_list);
}


void LatPacker::onSend(LatParticleData *buffer, const unsigned long send_len,
                       const int dimension, const int direction) {
    double offset[DIMENSION] = {0.0, 0.0, 0.0};
    // 当原子要跨越周期性边界, 原子坐标必须要做出调整
    setOffset(offset, dimension, direction);
    const int index = 2 * dimension + direction;
    pack::pack_send(dimension, send_len, offset, atom_list, buffer, send_list[index]);
}

void LatPacker::onReceive(LatParticleData *buffer, const unsigned long receive_len,
                          const int dimension, const int direction) {
    pack::unpack_recv(dimension, direction, receive_len, atom_list, buffer, receive_list);
}
