//
// Created by genshen on 2019-03-17.
//

#include <vector>
#include <logs/logs.h>
#include "atom/atom_element.h"
#include "atom/atom_list.h"
#include "lat_particle_packer.h"

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
    std::vector<_type_atom_id> &sendlist = send_list[index];
    _type_atom_id local_id;
    for (int i = 0; i < send_len; i++) {
        local_id = sendlist[i];
        AtomElement &atom = atom_list.getAtomEleByLinearIndex(local_id);
        // for ghost atoms, we just care their position and atom type(EamParser calculating), so positions and types are enough.
        buffer[i].type = atom.type; // fixme --type
        buffer[i].r[0] = atom.x[0] + offset[0];
        buffer[i].r[1] = atom.x[1] + offset[1];
        buffer[i].r[2] = atom.x[2] + offset[2];
    }
}

void LatPackerFirst::onReceive(LatParticleData *buffer, const unsigned long receive_len,
                               const int dimension, const int direction) {
    //将收到的粒子位置信息加到对应存储位置上
    const int d = dimension;
    const int n = receive_len;
    const _type_lattice_size (&ghost)[DIMENSION] = domain.dbx_lattice_size_ghost;
    const _type_lattice_size (&box)[DIMENSION] = domain.dbx_lattice_size_sub_box;
    const _type_lattice_size (&ext)[DIMENSION] = domain.dbx_lattice_size_ghost_extended;
    LatParticleData *buf = buffer;
    std::vector<std::vector<_type_atom_id> > &recvlist = receive_list;

    int xstart, ystart, zstart;
    int xstop, ystop, zstop;
    int recv_idnex = 2 * d + direction;
    int m = 0;
    if (d == 0) {
        if (direction == 0) { // mirror with send.
            xstart = ghost[0] + box[0];
            xstop = ext[0];
        } else {
            xstart = 0;
            xstop = ghost[0];
        }
        ystart = ghost[1];
        ystop = ystart + box[1];
        zstart = ghost[2];
        zstop = zstart + box[2];
        for (int k = zstart; k < zstop; k++) {
            for (int j = ystart; j < ystop; j++) {
                for (int i = xstart; i < xstop; i++) {
//                        kk = IndexOf3DIndex(i, j, k);
                    AtomElement &atom_ = atom_list.getAtomEleByGhostIndex(i, j, k);
                    atom_.type = buf[m].type;
                    atom_.x[0] = buf[m].r[0];
                    atom_.x[1] = buf[m].r[1];
                    atom_.x[2] = buf[m++].r[2];
                    recvlist[recv_idnex].push_back(atom_list.IndexOf3DIndex(i, j, k));
                }
            }
        }
        if (n != recvlist[recv_idnex].size()) { // todo error handling in dataReuse feature.
            kiwi::logs::e("unpack_recvfirst", "received data size does not match the Mpi_Proble size.\n");
        }

    } else if (d == 1) {
        if (direction == 0) {
            ystart = ghost[1] + box[1];
            ystop = ext[1];
        } else {
            ystart = 0;
            ystop = ghost[1];
        }
        xstart = 0;
        xstop = ext[0];
        zstart = ghost[2];
        zstop = zstart + box[2];

        for (int k = zstart; k < zstop; k++) {
            for (int j = ystart; j < ystop; j++) {
                for (int i = xstart; i < xstop; i++) {
                    AtomElement &atom_ = atom_list.getAtomEleByGhostIndex(i, j, k);
                    atom_.type = buf[m].type;
                    atom_.x[0] = buf[m].r[0];
                    atom_.x[1] = buf[m].r[1];
                    atom_.x[2] = buf[m++].r[2];
                    recvlist[recv_idnex].push_back(atom_list.IndexOf3DIndex(i, j, k));
                }
            }
        }
        if (n != recvlist[recv_idnex].size()) { // todo error handling in dataReuse feature.
            kiwi::logs::e("unpack_recvfirst", "received data size does not match the Mpi_Proble size.\n");
        }
    } else {
        if (direction == 0) {
            zstart = ghost[2] + box[2];
            zstop = ext[2];
        } else {
            zstart = 0;
            zstop = ghost[2];
        }
        xstart = 0;
        xstop = ext[0];
        ystart = 0;
        ystop = ext[1];
        for (int k = zstart; k < zstop; k++) {
            for (int j = ystart; j < ystop; j++) {
                for (int i = xstart; i < xstop; i++) {
//                        kk = IndexOf3DIndex(i, j, k);
                    AtomElement &atom_ = atom_list.getAtomEleByGhostIndex(i, j, k);
                    atom_.type = buf[m].type;
                    atom_.x[0] = buf[m].r[0];
                    atom_.x[1] = buf[m].r[1];
                    atom_.x[2] = buf[m++].r[2];
                    recvlist[recv_idnex].push_back(atom_list.IndexOf3DIndex(i, j, k));
                }
            }
        }
        if (n != recvlist[recv_idnex].size()) { // todo error handling in dataReuse feature.
            kiwi::logs::e("unpack_recvfirst", "received data size does not match the Mpi_Proble size.\n");
        }
    }
}


void LatPacker::onSend(LatParticleData *buffer, const unsigned long send_len,
                       const int dimension, const int direction) {
    double offset[DIMENSION] = {0.0, 0.0, 0.0};
    // 当原子要跨越周期性边界, 原子坐标必须要做出调整
    setOffset(offset, dimension, direction);
    const int index = 2 * dimension + direction;
    std::vector<_type_atom_id> &sendlist = send_list[index];
    _type_atom_id local_id;
    for (int i = 0; i < send_len; i++) {
        local_id = sendlist[i];
        AtomElement &atom = atom_list.getAtomEleByLinearIndex(local_id);
        // for ghost atoms, we just care their position and atom type(EamParser calculating), so positions and types are enough.
        buffer[i].type = atom.type; // fixme --type
        buffer[i].r[0] = atom.x[0] + offset[0];
        buffer[i].r[1] = atom.x[1] + offset[1];
        buffer[i].r[2] = atom.x[2] + offset[2];
    }
}

void LatPacker::onReceive(LatParticleData *buffer, const unsigned long receive_len,
                          const int dimension, const int direction) {
    int const n = receive_len;
    LatParticleData *buf = buffer;
    std::vector<std::vector<_type_atom_id>> &recvlist = receive_list;

    const int list_index = 2 * dimension + direction; // fixme Flip the direction

    long kk;
    for (int i = 0; i < n; i++) {
        kk = recvlist[list_index][i];
        AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(kk);
        atom_.type = buf[i].type;
        atom_.x[0] = buf[i].r[0];
        atom_.x[1] = buf[i].r[1];
        atom_.x[2] = buf[i].r[2];
    }
}
