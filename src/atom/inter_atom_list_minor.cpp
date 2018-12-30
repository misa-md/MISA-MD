//
// Created by genshen on 2018-12-28.
//

#include "inter_atom_list.h"
#include "pack/particledata.h"
#include "../domain/domain.h"
#include "../utils/mpi_domain.h"
#include "../utils/mpi_data_types.h"

void InterAtomList::exchangeInter(Domain *p_domain) {
    // 发送、接收数据缓冲区
    int numPartsToSend[DIMENSION][2];
    int numPartsToRecv[DIMENSION][2];
    particledata *sendbuf[2];
    particledata *recvbuf[2];

    // MPI通信状态和请求
    MPI_Status status;
    MPI_Status send_statuses[DIMENSION][2];
    MPI_Status recv_statuses[DIMENSION][2];
    MPI_Request send_requests[DIMENSION][2];
    MPI_Request recv_requests[DIMENSION][2];

    int direction;
    // 找到要发送给邻居的原子
    for (unsigned short d = 0; d < DIMENSION; d++) {
        for (direction = LOWER; direction <= HIGHER; direction++) {
            numPartsToSend[d][direction] = countExSendNum(p_domain, d, direction);
            sendbuf[direction] = new particledata[numPartsToSend[d][direction]];
            packExInterToSend(p_domain, sendbuf[direction], d, direction);
        }

        // 与上下邻居通信
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numsend = numPartsToSend[d][direction];
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

unsigned long InterAtomList::countExSendNum(Domain *p_domain, int dimension, int direction) {
    unsigned long i = 0;
    if (direction == 0) {
        for (AtomElement &inter_ref :inter_list) {
            // we assume that, a atom cannot cross 2 or more than 2 sub-boxes
            if (inter_ref.x[dimension] < p_domain->meas_sub_box_region.low[dimension]) {
                i++;
            }
        }
    } else {
        for (AtomElement &inter_ref :inter_list) {
            if (inter_ref.x[dimension] >= p_domain->meas_sub_box_region.high[dimension]) {
                i++;
            }
        }
    }
    return i;
}

void InterAtomList::packExInterToSend(Domain *p_domain, particledata *buf, int dimension, int direction) {
    if (direction == 0) {
        unsigned long i = 0; // todo type
        for (_type_inter_list::iterator inter_it = inter_list.begin();
             inter_it != inter_list.end();) {
            // we assume that, a atom cannot cross 2 or more than 2 sub-boxes
            if ((*inter_it).x[dimension] < p_domain->meas_sub_box_region.low[dimension]) {
                buf[i].id = inter_it->id;
                buf[i].type = inter_it->type;
                buf[i].r[0] = inter_it->x[0];
                buf[i].r[1] = inter_it->x[1];
                buf[i].r[2] = inter_it->x[2];
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
    } else {
        unsigned long i = 0;
        for (_type_inter_list::iterator inter_it = inter_list.begin();
             inter_it != inter_list.end();) {
            // we assume that, a atom cannot cross 2 or more than 2 sub-boxes
            if (inter_it->x[dimension] >= p_domain->meas_sub_box_region.high[dimension]) {
                buf[i].id = inter_it->id;
                buf[i].type = inter_it->type;
                buf[i].r[0] = inter_it->x[0];
                buf[i].r[1] = inter_it->x[1];
                buf[i].r[2] = inter_it->x[2];
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
}

void InterAtomList::unpackExInterRecv(int d, int n, const double *lower, const double *upper, particledata *buf) {
    AtomElement atom;
    atom_type::atom_type type;
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
