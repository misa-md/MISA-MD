//
// Created by genshen on 2018-12-28.
//

#include <logs/logs.h>
#include "inter_atom_list.h"
#include "pack/particledata.h"
#include "../domain/domain.h"
#include "../utils/mpi_domain.h"
#include "../utils/mpi_data_types.h"

void InterAtomList::exchangeInter(Domain *p_domain) {
    // 发送、接收数据缓冲区
    int num_parts_to_send[DIMENSION][2];
    int numPartsToRecv[DIMENSION][2];
    particledata *sendbuf[2];
    particledata *recvbuf[2];

    // MPI通信状态和请求
    MPI_Status status;
    MPI_Status send_statuses[DIMENSION][2];
    MPI_Status recv_statuses[DIMENSION][2];
    MPI_Request send_requests[DIMENSION][2];
    MPI_Request recv_requests[DIMENSION][2];

    // pre calculate the number of send data
    countExSendNum(p_domain, num_parts_to_send);

    box::_type_flag_32 out_box_flags[DIMENSION][2] = {
            {box::OUT_BOX_X_LITTER, box::OUT_BOX_X_BIG},
            {box::OUT_BOX_Y_LITTER, box::OUT_BOX_Y_BIG},
            {box::OUT_BOX_Z_LITTER, box::OUT_BOX_Z_BIG},
    };
    int direction;
    double offset[DIMENSION] = {0.0};
    // 找到要发送给邻居的原子
    for (unsigned short d = 0; d < DIMENSION; d++) {
        for (direction = LOWER; direction <= HIGHER; direction++) {
            for (double &_d : offset) {
                _d = 0.0;
            }
            // periodic boundary
            if (p_domain->grid_coord_sub_box[d] == 0 && direction == LOWER) {
                offset[d] = p_domain->meas_global_length[d];
            }
            if (p_domain->grid_coord_sub_box[d] == p_domain->grid_size[d] - 1 && direction == HIGHER) {
                offset[d] = -((p_domain->meas_global_length[d]));
            }

            sendbuf[direction] = new particledata[num_parts_to_send[d][direction]];
            packExInterToSend(p_domain, sendbuf[direction], d, direction, out_box_flags, offset);
        }

        // 与上下邻居通信
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numsend = num_parts_to_send[d][direction];
            int numrecv;

            MPI_Isend(sendbuf[direction], numsend, mpi_types::_mpi_Particle_data,
                      p_domain->rank_id_neighbours[d][direction], 99,
                      MPIDomain::sim_processor.comm, &send_requests[d][direction]);
            MPI_Probe(p_domain->rank_id_neighbours[d][(direction + 1) % 2], 99, MPIDomain::sim_processor.comm,
                      &status);//测试邻居是否有信息发送给本地
            MPI_Get_count(&status, mpi_types::_mpi_Particle_data, &numrecv);//得到要接收的粒子数目
            // 初始化接收缓冲区
            //依据得到发送方要发送粒子信息大小，初始化接收缓冲区
            recvbuf[direction] = new particledata[numrecv];
            numPartsToRecv[d][direction] = numrecv;
            MPI_Irecv(recvbuf[direction], numrecv, mpi_types::_mpi_Particle_data,
                      p_domain->rank_id_neighbours[d][(direction + 1) % 2],
                      99, MPIDomain::sim_processor.comm, &recv_requests[d][direction]);
        } // todo combine this two loops.

        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numrecv = numPartsToRecv[d][direction];
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);
            //将收到的粒子位置信息加到对应存储位置上
            unpackExInterRecv(d, numrecv,
                              p_domain->meas_sub_box_region.low,
                              p_domain->meas_sub_box_region.high,
                              recvbuf[direction]);
//            _atom->unpackExInterRecv(d, numrecv, recvbuf[direction]);

            // 释放buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}

void InterAtomList::countExSendNum(Domain *p_domain, int n_to_send[DIMENSION][2]) {
    box::_type_flag_32 flag;
    // clear data.
    for (int d = 0; d < DIMENSION; d++) {
        for (int dir = LOWER; dir <= HIGHER; dir++) {
            n_to_send[d][dir] = 0;
        }
    }

    for (AtomElement &inter_ref :inter_list) {
        // we assume that, a atom cannot cross 2 or more than 2 sub-boxes
        flag = isOutBox(inter_ref, p_domain);
        if (flag & box::OUT_BOX_X_LITTER) {
            n_to_send[0][LOWER]++;
        } else if (flag & box::OUT_BOX_X_BIG) {
            n_to_send[0][HIGHER]++;
        } else if (flag & box::OUT_BOX_Y_LITTER) {
            n_to_send[1][LOWER]++;
        } else if (flag & box::OUT_BOX_Y_BIG) {
            n_to_send[1][HIGHER]++;
        } else if (flag & box::OUT_BOX_Z_LITTER) {
            n_to_send[2][LOWER]++;
        } else if (flag & box::OUT_BOX_Z_BIG) {
            n_to_send[2][HIGHER]++;
        }
    }
}

void InterAtomList::packExInterToSend(Domain *p_domain, particledata *buf, int dimension, int direction,
                                      box::_type_flag_32 excepted_flag[DIMENSION][2], double offset[DIMENSION]) {
    unsigned long i = 0; // todo type
    for (_type_inter_list::iterator inter_it = inter_list.begin(); inter_it != inter_list.end();) {
        // we assume that, a atom cannot cross 2 or more than 2 sub-boxes
        if (isOutBox(*inter_it, p_domain) & excepted_flag[dimension][direction]) {
            buf[i].id = inter_it->id;
            buf[i].type = inter_it->type;
            buf[i].r[0] = inter_it->x[0] + offset[0];
            buf[i].r[1] = inter_it->x[1] + offset[1];
            buf[i].r[2] = inter_it->x[2] + offset[2];
            buf[i].v[0] = inter_it->v[0];
            buf[i].v[1] = inter_it->v[1];
            buf[i].v[2] = inter_it->v[2];
            // remove the inter atom.
            // exchange the atom at end of vector to the position of atom (inter[j]) te be removed .
            inter_it = inter_list.erase(inter_it);
            nlocalinter--;
            i++;
        } else {
            inter_it++;
        }
    }
}

void InterAtomList::unpackExInterRecv(int d, int n, const double *lower, const double *upper, particledata *buf) {
    AtomElement atom;
    for (int i = 0; i < n; i++) {
        atom.id = buf[i].id;
        atom.type = buf[i].type;
        atom.x[0] = buf[i].r[0];
        atom.x[1] = buf[i].r[1];
        atom.x[2] = buf[i].r[2];
        atom.v[0] = buf[i].v[0];
        atom.v[1] = buf[i].v[1];
        atom.v[2] = buf[i].v[2];
        if (atom.x[d] >= lower[d] && atom.x[d] < upper[d]) {
            inter_list.push_back(atom);
            nlocalinter++;
        } else {
            // todo waring
        }
    }
}
