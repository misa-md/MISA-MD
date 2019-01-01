//
// Created by genshen on 5/19/18.
//

#include <cmath>
#include "inter_atom_list.h"
#include "../domain/domain.h"
#include "../utils/mpi_domain.h"
#include "../utils/mpi_data_types.h"

InterAtomList::InterAtomList() : nlocalinter(0), nghostinter(0),
                                 intersendlist(6), interrecvlist(6) {}

void InterAtomList::appendInter(_type_atom_id atom_id) {

}

void InterAtomList::addInterAtom(AtomElement &atom) {
    inter_list.push_back(atom);
    nlocalinter++;
    // todo set df,f,rhointer to 0.
}

const box::_type_flag_32 InterAtomList::isOutBox(const AtomElement &src_atom, const Domain *p_domain) {
    auto j = static_cast<_type_atom_index>(lround(src_atom.x[0] * 2 / p_domain->lattice_const));
    auto k = static_cast<_type_atom_index>(lround(src_atom.x[1] * 2 / p_domain->lattice_const));
    auto l = static_cast<_type_atom_index>(lround(src_atom.x[2] * 2 / p_domain->lattice_const));
    k = k < 0 ? k / 2 - 1 : k / 2; // for example: right k=-1, then -1/2 = 0, but we expect left k=-1
    l = l < 0 ? l / 2 - 1 : l / 2;
    j -= 2 * p_domain->lattice_coord_sub_box_region.x_low;
    k -= p_domain->lattice_coord_sub_box_region.y_low;
    l -= p_domain->lattice_coord_sub_box_region.z_low;

    box::_type_flag_32 flag = box::IN_BOX;
    if (j < 0) {
        flag |= box::OUT_BOX_X_LITTER;
    } else if (j >= 2 * p_domain->lattice_size_sub_box[0]) {
        flag |= box::OUT_BOX_X_BIG;
    }
    if (k < 0) {
        flag |= box::OUT_BOX_Y_LITTER;
    } else if (k >= p_domain->lattice_size_sub_box[1]) {
        flag |= box::OUT_BOX_Y_BIG;
    }
    if (l < 0) {
        flag |= box::OUT_BOX_Z_LITTER;
    } else if (l >= p_domain->lattice_size_sub_box[2]) {
        flag |= box::OUT_BOX_Z_BIG;
    }
    return flag;
}

void InterAtomList::borderInter(Domain *p_domain) {
    intersendlist.clear();
    interrecvlist.clear();
    intersendlist.resize(6);
    interrecvlist.resize(6);
    // 发送、接收数据缓冲区
    int numPartsToSend[DIMENSION][2];
    int numPartsToRecv[DIMENSION][2];
    LatParticleData *sendbuf[2];
    LatParticleData *recvbuf[2];

    // MPI通信状态和请求
    MPI_Status status;
    MPI_Status send_statuses[DIMENSION][2];
    MPI_Status recv_statuses[DIMENSION][2];
    MPI_Request send_requests[DIMENSION][2];
    MPI_Request recv_requests[DIMENSION][2];

    int direction;
    int iswap = 0;
    int jswap = 0;
    for (unsigned short d = 0; d < DIMENSION; d++) {
        double offsetLower[DIMENSION];
        double offsetHigher[DIMENSION];
        offsetLower[d] = 0.0;
        offsetHigher[d] = 0.0;

        // 进程在左侧边界
        if (p_domain->grid_coord_sub_box[d] == 0) {
            offsetLower[d] = p_domain->meas_global_length[d];
        }
        // 进程在右侧边界
        if (p_domain->grid_coord_sub_box[d] == p_domain->grid_size[d] - 1) {
            offsetHigher[d] = -((p_domain->meas_global_length[d]));
        }

        for (direction = LOWER; direction <= HIGHER; direction++) {
            // 找到要发送给邻居的原子
            getIntertosend(p_domain, d, direction, p_domain->meas_ghost_length[d], intersendlist[iswap]);

            double shift = 0.0;
            if (direction == LOWER) {
                shift = offsetLower[d];
            }
            if (direction == HIGHER) {
                shift = offsetHigher[d];
            }

            // 初始化发送缓冲区
            numPartsToSend[d][direction] = intersendlist[iswap].size();
            sendbuf[direction] = new LatParticleData[numPartsToSend[d][direction]];
            pack_bordersend(d, numPartsToSend[d][direction],
                            intersendlist[iswap++], sendbuf[direction], shift);
//            _atom->pack_bordersend(d, numPartsToSend[d][direction], intersendlist[iswap++], sendbuf[direction], shift);
        }

        // 与上下邻居通信
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numsend = numPartsToSend[d][direction];
            int numrecv;

            MPI_Isend(sendbuf[direction], numsend,
                      mpi_types::_mpi_latParticle_data,
                      p_domain->rank_id_neighbours[d][direction],
                      99, MPIDomain::sim_processor.comm, &send_requests[d][direction]);
            MPI_Probe(p_domain->rank_id_neighbours[d][(direction + 1) % 2], 99, MPIDomain::sim_processor.comm,
                      &status);//测试邻居是否有信息发送给本地
            MPI_Get_count(&status, mpi_types::_mpi_latParticle_data, &numrecv);//得到要接收的粒子数目
            // 初始化接收缓冲区
            //依据得到发送方要发送粒子信息大小，初始化接收缓冲区
            recvbuf[direction] = new LatParticleData[numrecv];
            numPartsToRecv[d][direction] = numrecv;
            MPI_Irecv(recvbuf[direction], numrecv, mpi_types::_mpi_latParticle_data,
                      p_domain->rank_id_neighbours[d][(direction + 1) % 2], 99,
                      MPIDomain::sim_processor.comm, &recv_requests[d][direction]);
        }

        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numrecv = numPartsToRecv[d][direction];
            interrecvlist[jswap].resize(numrecv);
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);

            //将收到的粒子位置信息加到对应存储位置上
            unpack_borderrecv(numrecv,
                              p_domain->meas_ghost_region.low,
                              p_domain->meas_ghost_region.high,
                              recvbuf[direction], interrecvlist[jswap++]);
//            _atom->unpack_borderrecv(numrecv, recvbuf[direction], interrecvlist[jswap++]);
            // 释放buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}

void InterAtomList::pack_bordersend(int dimension, int n, std::vector<AtomElement *> &sendlist,
                                    LatParticleData *buf, double shift) {
    // fixme pack atom id?
    switch (dimension) {
        case 0:
            for (int i = 0; i < n; i++) {
                buf[i].type = sendlist[i]->type;
                buf[i].r[0] = sendlist[i]->x[0] + shift;
                buf[i].r[1] = sendlist[i]->x[1];
                buf[i].r[2] = sendlist[i]->x[2];
            }
            break;
        case 1:
            for (int i = 0; i < n; i++) {
                buf[i].type = sendlist[i]->type;
                buf[i].r[0] = sendlist[i]->x[0];
                buf[i].r[1] = sendlist[i]->x[1] + shift;
                buf[i].r[2] = sendlist[i]->x[2];
            }
            break;
        case 2:
            for (int i = 0; i < n; i++) {
                buf[i].type = sendlist[i]->type;
                buf[i].r[0] = sendlist[i]->x[0];
                buf[i].r[1] = sendlist[i]->x[1];
                buf[i].r[2] = sendlist[i]->x[2] + shift;
            }
            break;
        default:
            break;
    }
}

void InterAtomList::unpack_borderrecv(int n, const double lower[DIMENSION], const double upper[DIMENSION],
                                      LatParticleData *buf, std::vector<AtomElement *> &recvlist) {
    atom_type::atom_type type;
    // std::vector<double> xtemp(3);
    AtomElement ele;
    for (int i = 0; i < n; i++) {
        ele.type = buf[i].type;
        ele.x[0] = buf[i].r[0];
        ele.x[1] = buf[i].r[1];
        ele.x[2] = buf[i].r[2];
        if (ele.x[0] >= lower[0] && ele.x[0] < upper[0] &&
            ele.x[1] >= lower[1] && ele.x[1] < upper[1] &&
            ele.x[2] >= lower[2] && ele.x[2] < upper[2]) {
            inter_ghost_list.push_back(ele);
            nghostinter++;
        } else {
            // todo warning
            recvlist[i] = nullptr;
        }
    }
}

void InterAtomList::getIntertosend(Domain *p_domain, int d, int direction, double ghostlengh,
                                   std::vector<AtomElement *> &sendlist) {
    double low, high;
    if (d == 0) {
        if (direction == 0) {
            low = p_domain->meas_sub_box_region.x_low;
            high = p_domain->meas_sub_box_region.x_low + ghostlengh;
            for (AtomElement &inter_ref : inter_list) {
                if (inter_ref.x[0] < high && inter_ref.x[0] >= low) {
                    sendlist.push_back(&inter_ref);
                }
            }
            for (AtomElement &ghost_ref :inter_ghost_list) {
                if (ghost_ref.x[0] < high && ghost_ref.x[0] >= low) {
                    sendlist.push_back(&ghost_ref);
                }
            }
        } else {
            low = p_domain->meas_sub_box_region.x_high - ghostlengh;
            high = p_domain->meas_sub_box_region.x_high;
            for (AtomElement &inter_ref :inter_list) {
                if (inter_ref.x[0] <= high && inter_ref.x[0] > low) {
                    sendlist.push_back(&inter_ref);
                }
            }
            for (AtomElement &ghost_ref :inter_ghost_list) {
                if (ghost_ref.x[0] <= high && ghost_ref.x[0] > low) {
                    sendlist.push_back(&ghost_ref);
                }
            }
        }
    } else if (d == 1) {
        if (direction == 0) {
            low = p_domain->meas_sub_box_region.y_low;
            high = p_domain->meas_sub_box_region.y_low + ghostlengh;
            for (AtomElement &inter_ref :inter_list) {
                if (inter_ref.x[1] < high && inter_ref.x[1] >= low) {
                    sendlist.push_back(&inter_ref);
                }
            }
            for (AtomElement &ghost_ref :inter_ghost_list) {
                if (ghost_ref.x[1] < high && ghost_ref.x[1] >= low) {
                    sendlist.push_back(&ghost_ref);
                }
            }
        } else {
            low = p_domain->meas_sub_box_region.y_high - ghostlengh;
            high = p_domain->meas_sub_box_region.y_high;
            for (AtomElement &inter_ref :inter_list) {
                if (inter_ref.x[1] <= high && inter_ref.x[1] > low) {
                    sendlist.push_back(&inter_ref);
                }
            }
            for (AtomElement &ghost_ref :inter_ghost_list) {
                if (ghost_ref.x[1] <= high && ghost_ref.x[1] > low) {
                    sendlist.push_back(&ghost_ref);
                }
            }
        }
    } else {
        if (direction == 0) {
            low = p_domain->meas_sub_box_region.z_low;
            high = p_domain->meas_sub_box_region.z_low + ghostlengh;
            for (AtomElement &inter_ref :inter_list) {
                if (inter_ref.x[2] < high && inter_ref.x[2] >= low) {
                    sendlist.push_back(&inter_ref);
                }
            }
            for (AtomElement &ghost_ref :inter_ghost_list) {
                if (ghost_ref.x[2] < high && ghost_ref.x[2] >= low) {
                    sendlist.push_back(&ghost_ref);
                }
            }
        } else {
            low = p_domain->meas_sub_box_region.z_high - ghostlengh;
            high = p_domain->meas_sub_box_region.z_high;
            for (AtomElement &inter_ref :inter_list) {
                if (inter_ref.x[2] <= high && inter_ref.x[2] > low) {
                    sendlist.push_back(&inter_ref);
                }
            }
            for (AtomElement &ghost_ref :inter_ghost_list) {
                if (ghost_ref.x[2] <= high && ghost_ref.x[2] > low) {
                    sendlist.push_back(&ghost_ref);
                }
            }
        }
    }
}
