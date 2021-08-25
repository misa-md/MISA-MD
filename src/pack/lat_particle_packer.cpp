//
// Created by genshen on 2019-03-17.
//

#include <vector>
#include <logs/logs.h>
#include "atom/atom_element.h"
#include "atom/atom_list.h"
#include "lat_particle_packer.h"

LatParticlePacker::LatParticlePacker(const comm::BccDomain &domain, AtomList &atom_list,
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
    if (domain.grid_coord[dimension] == 0 && direction == comm::DIR_LOWER) {
        offset[dimension] = domain.meas_global_length[dimension];
    }
    // 进程在右侧边界
    if (domain.grid_coord[dimension] == domain.grid_size[dimension] - 1 && direction == comm::DIR_HIGHER) {
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
        MD_LOAD_ATOM_VAR(atom, (&atom_list), local_id);

        // for ghost atoms, we just care their position and atom type(EamParser calculating), so positions and types are enough.
        buffer[i].type = MD_GET_ATOM_TYPE(atom, local_id); // fixme --type
#ifdef DEV_MD_COMM_INC_ATOM_ID
        buffer[i].id = MD_GET_ATOM_ID(atom, local_id);
#endif
        buffer[i].r[0] = MD_GET_ATOM_X(atom, local_id, 0) + offset[0];
        buffer[i].r[1] = MD_GET_ATOM_X(atom, local_id, 1) + offset[1];
        buffer[i].r[2] = MD_GET_ATOM_X(atom, local_id, 2) + offset[2];
    }
}

void LatPackerFirst::onReceive(LatParticleData *buffer, const unsigned long receive_len,
                               const int dimension, const int direction) {
    //将收到的粒子位置信息加到对应存储位置上
    const _type_lattice_size (&ghost)[DIMENSION] = domain.dbx_lattice_size_ghost;
    const _type_lattice_size (&box)[DIMENSION] = domain.dbx_sub_box_lattice_size;
    const _type_lattice_size (&ext)[DIMENSION] = domain.dbx_ghost_extended_lattice_size;

    int xstart, ystart, zstart;
    int xstop, ystop, zstop;
    int recv_index = 2 * dimension + direction;
    int m = 0;
    if (dimension == 0) {
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
                    const _type_atom_index gid = atom_list._atoms.getAtomIndex(i, j, k);
                    MD_LOAD_ATOM_VAR(atom_, (&atom_list), gid);
                    MD_SET_ATOM_TYPE(atom_, gid, buffer[m].type);
#ifdef DEV_MD_COMM_INC_ATOM_ID
                    atom_.id = buffer[m].id;
#endif
                    MD_SET_ATOM_X(atom_, gid, 0, buffer[m].r[0]);
                    MD_SET_ATOM_X(atom_, gid, 1, buffer[m].r[1]);
                    MD_SET_ATOM_X(atom_, gid, 2, buffer[m++].r[2]);
                    receive_list[recv_index].push_back(atom_list.lattice.IndexOf3DIndex(i, j, k));
                }
            }
        }
        if (receive_len != receive_list[recv_index].size()) { // todo error handling in dataReuse feature.
            kiwi::logs::e("unpack_recvfirst",
                          "received data size does not match the MPI_Proble size，expected {}, but got {}.\n",
                          receive_len,
                          receive_list[recv_index].size());
        }

    } else if (dimension == 1) {
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
                    const _type_atom_index gid = atom_list._atoms.getAtomIndex(i, j, k);
                    MD_LOAD_ATOM_VAR(atom_, (&atom_list), gid);
                    MD_SET_ATOM_TYPE(atom_, gid, buffer[m].type);
#ifdef DEV_MD_COMM_INC_ATOM_ID
                    atom_.id = buffer[m].id;
#endif
                    MD_SET_ATOM_X(atom_, gid, 0, buffer[m].r[0]);
                    MD_SET_ATOM_X(atom_, gid, 1, buffer[m].r[1]);
                    MD_SET_ATOM_X(atom_, gid, 2, buffer[m++].r[2]);
                    receive_list[recv_index].push_back(atom_list.lattice.IndexOf3DIndex(i, j, k));
                }
            }
        }
        if (receive_len != receive_list[recv_index].size()) { // todo error handling in dataReuse feature.
            kiwi::logs::e("unpack_recvfirst",
                          "received data size does not match the MPI_Proble size，expected {}, but got {}.\n",
                          receive_len,
                          receive_list[recv_index].size());
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
                    const _type_atom_index gid = atom_list._atoms.getAtomIndex(i, j, k);
                    MD_LOAD_ATOM_VAR(atom_, (&atom_list), gid);
                    MD_SET_ATOM_TYPE(atom_, gid, buffer[m].type);
#ifdef DEV_MD_COMM_INC_ATOM_ID
                    atom_.id = buffer[m].id;
#endif
                    MD_SET_ATOM_X(atom_, gid, 0, buffer[m].r[0]);
                    MD_SET_ATOM_X(atom_, gid, 1, buffer[m].r[1]);
                    MD_SET_ATOM_X(atom_, gid, 2, buffer[m++].r[2]);
                    receive_list[recv_index].push_back(atom_list.lattice.IndexOf3DIndex(i, j, k));
                }
            }
        }
        if (receive_len != receive_list[recv_index].size()) { // todo error handling in dataReuse feature.
            kiwi::logs::e("unpack_recvfirst",
                          "received data size does not match the MPI_Proble size，expected {}, but got {}.\n",
                          receive_len,
                          receive_list[recv_index].size());
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
        MD_LOAD_ATOM_VAR(atom, (&atom_list), local_id);
        // for ghost atoms, we just care their position and atom type(EamParser calculating), so positions and types are enough.
        buffer[i].type = MD_GET_ATOM_TYPE(atom, local_id); // fixme --type
#ifdef DEV_MD_COMM_INC_ATOM_ID
        buffer[i].id = atom.id;
#endif
        buffer[i].r[0] = MD_GET_ATOM_X(atom, local_id, 0) + offset[0];
        buffer[i].r[1] = MD_GET_ATOM_X(atom, local_id, 1) + offset[1];
        buffer[i].r[2] = MD_GET_ATOM_X(atom, local_id, 2) + offset[2];
    }
}

void LatPacker::onReceive(LatParticleData *buffer, const unsigned long receive_len,
                          const int dimension, const int direction) {
    const int list_index = 2 * dimension + direction; // fixme Flip the direction
    long kk;
    for (int i = 0; i < receive_len; i++) {
        kk = receive_list[list_index][i];
        MD_LOAD_ATOM_VAR(atom_, (&atom_list), kk);

        MD_SET_ATOM_TYPE(atom_, kk, buffer[i].type);
#ifdef DEV_MD_COMM_INC_ATOM_ID
        atom_.id = buffer[i].id;
#endif
        MD_SET_ATOM_X(atom_, kk, 0, buffer[i].r[0]);
        MD_SET_ATOM_X(atom_, kk, 1, buffer[i].r[1]);
        MD_SET_ATOM_X(atom_, kk, 2, buffer[i].r[2]);
    }
}
