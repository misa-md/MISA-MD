//
// Created by genshen on 5/19/18.
//

#include "inter_atom_list.h"
#include "utils/mpi_domain.h"
#include "utils/mpi_data_types.h"
#include "../domain.h"

InterAtomList::InterAtomList() : nlocalinter(0), nghostinter(0) {

}

void InterAtomList::appendInter(_type_atom_id atom_id) {

}

void InterAtomList::addInterAtom(AtomElement &atom) {
    inter_list.push_back(atom);
    nlocalinter++;
    // todo set df,f,rhointer to 0.
}


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
            numPartsToSend[d][direction] = getintersendnum(p_domain, d, direction);
            sendbuf[direction] = new particledata[numPartsToSend[d][direction]];
            pack_intersend(interbuf, sendbuf[direction]);
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
            unpack_interrecv(d, numrecv,
                             p_domain->meas_sub_box_lower_bounding,
                             p_domain->meas_sub_box_upper_bounding,
                             recvbuf[direction]);
//            _atom->unpack_interrecv(d, numrecv, recvbuf[direction]);

            // 释放buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
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
                              p_domain->meas_ghost_lower_bounding,
                              p_domain->meas_ghost_upper_bounding,
                              recvbuf[direction], interrecvlist[jswap++]);
//            _atom->unpack_borderrecv(numrecv, recvbuf[direction], interrecvlist[jswap++]);
            // 释放buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}

void InterAtomList::pack_intersend(std::vector<unsigned long> interbuf, particledata *buf) {
    int j;
    for (int i = 0; i < interbuf.size(); i++) {
        j = interbuf[i];
        buf[i].id = inter->idinter[j];
        buf[i].type = inter->typeinter[j];
        buf[i].r[0] = inter->xinter[j][0];
        buf[i].r[1] = inter->xinter[j][1];
        buf[i].r[2] = inter->xinter[j][2];
        buf[i].v[0] = inter->vinter[j][0];
        buf[i].v[1] = inter->vinter[j][1];
        buf[i].v[2] = inter->vinter[j][2];
        // remove the inter atom.
        // exchange the atom at end of vector to the position of atom (inter[j]) te be removed .
        inter->idinter[j] = inter->idinter[nlocalinter - 1];
        inter->typeinter[j] = inter->typeinter[nlocalinter - 1];
        inter->xinter[j][0] = inter->xinter[nlocalinter - 1][0];
        inter->xinter[j][1] = inter->xinter[nlocalinter - 1][1];
        inter->xinter[j][2] = inter->xinter[nlocalinter - 1][2];
        inter->vinter[j][0] = inter->vinter[nlocalinter - 1][0];
        inter->vinter[j][1] = inter->vinter[nlocalinter - 1][1];
        inter->vinter[j][2] = inter->vinter[nlocalinter - 1][2];
        nlocalinter--;
    }
}

void InterAtomList::unpack_interrecv(int d, int n,
                                     const double lower[DIMENSION], // p_domain->getMeasuredSubBoxLowerBounding(d)
                                     const double upper[DIMENSION], // p_domain->getMeasuredSubBoxUpperBounding(d)
                                     particledata *buf) {
    std::vector<double> xtemp(3);
    std::vector<double> vtemp(3);
    unsigned long id;
    atom_type::atom_type type;
    for (int i = 0; i < n; i++) {
        id = buf[i].id;
        type = buf[i].type;
        xtemp[0] = buf[i].r[0];
        xtemp[1] = buf[i].r[1];
        xtemp[2] = buf[i].r[2];
        vtemp[0] = buf[i].v[0];
        vtemp[1] = buf[i].v[1];
        vtemp[2] = buf[i].v[2];
        if (xtemp[d] >= lower[d] &&
            xtemp[d] < upper[d]) {
            if (nlocalinter == inter->xinter.size()) {
                inter->idinter.push_back(id);
                inter->typeinter.push_back(type);
                inter->xinter.push_back(xtemp);
                inter->vinter.push_back(vtemp);
                nlocalinter++;
                inter->finter.resize(nlocalinter, std::vector<double>(3));
                inter->rhointer.resize(nlocalinter);
                inter->dfinter.resize(nlocalinter);
            } else {
                if (inter->idinter.size() == nlocalinter) {
                    inter->idinter.push_back(id);
                } else {
                    inter->idinter[nlocalinter] = id;
                }
                inter->typeinter[nlocalinter] = type;
                inter->xinter[nlocalinter][0] = xtemp[0];
                inter->xinter[nlocalinter][1] = xtemp[1];
                inter->xinter[nlocalinter][2] = xtemp[2];
                if (nlocalinter == inter->vinter.size()) {
                    inter->vinter.push_back(vtemp);
                } else {
                    inter->vinter[nlocalinter][0] = vtemp[0];
                    inter->vinter[nlocalinter][1] = vtemp[1];
                    inter->vinter[nlocalinter][2] = vtemp[2];
                }
                nlocalinter++;
                inter->finter.resize(nlocalinter, std::vector<double>(3));
                inter->rhointer.resize(nlocalinter);
                inter->dfinter.resize(nlocalinter);
            }
        }
    }
}

void InterAtomList::pack_bordersend(int dimension, int n,
                                    std::vector<int> &sendlist, LatParticleData *buf, double shift) {
    int j;
    if (dimension == 0) {
        for (int i = 0; i < n; i++) {
            j = sendlist[i];
            buf[i].type = inter->typeinter[j];
            buf[i].r[0] = inter->xinter[j][0] + shift;
            buf[i].r[1] = inter->xinter[j][1];
            buf[i].r[2] = inter->xinter[j][2];
        }
    } else if (dimension == 1) {
        for (int i = 0; i < n; i++) {
            j = sendlist[i];
            buf[i].type = inter->typeinter[j];
            buf[i].r[0] = inter->xinter[j][0];
            buf[i].r[1] = inter->xinter[j][1] + shift;
            buf[i].r[2] = inter->xinter[j][2];
        }
    } else {
        for (int i = 0; i < n; i++) {
            j = sendlist[i];
            buf[i].type = inter->typeinter[j];
            buf[i].r[0] = inter->xinter[j][0];
            buf[i].r[1] = inter->xinter[j][1];
            buf[i].r[2] = inter->xinter[j][2] + shift;
        }
    }
}

void InterAtomList::unpack_borderrecv(int n,
                                      const double lower[DIMENSION], // p_domain->getMeasuredGhostLowerBounding(d)
                                      const double upper[DIMENSION], // p_domain->getMeasuredGhostUpperBounding(d)
                                      LatParticleData *buf, std::vector<int> &recvlist) {
    atom_type::atom_type type;
    std::vector<double> xtemp(3);
    for (int i = 0; i < n; i++) {
        type = buf[i].type;
        xtemp[0] = buf[i].r[0];
        xtemp[1] = buf[i].r[1];
        xtemp[2] = buf[i].r[2];
        if (xtemp[0] >= lower[0] &&
            xtemp[0] < upper[0] &&
            xtemp[1] >= lower[1] &&
            xtemp[1] < upper[1] &&
            xtemp[2] >= lower[2] &&
            xtemp[2] < upper[2]) {
            if (inter->xinter.size() == nlocalinter + nghostinter) {
                inter->typeinter.push_back(type);
                inter->xinter.push_back(xtemp);
                nghostinter++;
                recvlist[i] = nlocalinter + nghostinter - 1;
                inter->finter.resize(nlocalinter + nghostinter, std::vector<double>(3));
                inter->rhointer.resize(nlocalinter + nghostinter);
                inter->dfinter.resize(nlocalinter + nghostinter);
            } else {
                inter->typeinter[nlocalinter + nghostinter] = type;
                inter->xinter[nlocalinter + nghostinter][0] = xtemp[0];
                inter->xinter[nlocalinter + nghostinter][1] = xtemp[1];
                inter->xinter[nlocalinter + nghostinter][2] = xtemp[2];
                nghostinter++;
                recvlist[i] = nlocalinter + nghostinter - 1;
                inter->finter.resize(nlocalinter + nghostinter, std::vector<double>(3));
                inter->rhointer.resize(nlocalinter + nghostinter);
                inter->dfinter.resize(nlocalinter + nghostinter);
            }
        } else {
            recvlist[i] = -1;
        }
    }
}

unsigned long InterAtomList::getinteridsendsize() {
    return interbuf.size();
}

void InterAtomList::getIntertosend(Domain *p_domain, int d, int direction, double ghostlengh,
                                   std::vector<int> &sendlist) {
    double low, high;
    unsigned long i = 0;
    if (d == 0) {
        if (direction == 0) {
            low = p_domain->meas_sub_box_lower_bounding[0];
            high = p_domain->meas_sub_box_lower_bounding[0] + ghostlengh;
            // fixme fixme checkout ghost inters
            i = 0;
            for (AtomElement &inter_ref : inter_list) {
                if (inter_ref.x[0] < high && inter_ref.x[0] >= low) {
                    sendlist.push_back(i);
                }
                i++;
            }
            i = nLocalInter();
            for (AtomElement &ghost_ref :inter_ghost_list) {
                if (ghost_ref.x[0] < high && ghost_ref.x[0] >= low) {
                    sendlist.push_back(i);
                }
                i++;
            }
        } else {
            low = p_domain->meas_sub_box_upper_bounding[0] - ghostlengh;
            high = p_domain->meas_sub_box_upper_bounding[0];
            i = 0;
            for (AtomElement &inter_ref :inter_list) {
                if (inter_ref.x[0] <= high && inter_ref.x[0] > low) {
                    sendlist.push_back(i);
                }
                i++;
            }
            i = nLocalInter();
            for (AtomElement &ghost_ref :inter_ghost_list) {
                if (ghost_ref.x[0] <= high && ghost_ref.x[0] > low) {
                    sendlist.push_back(i);
                }
                i++;
            }
        }
    } else if (d == 1) {
        if (direction == 0) {
            low = p_domain->meas_sub_box_lower_bounding[1];
            high = p_domain->meas_sub_box_lower_bounding[1] + ghostlengh;
            i = 0;
            for (AtomElement &inter_ref :inter_list) {
                if (inter_ref.x[1] < high && inter_ref.x[1] >= low) {
                    sendlist.push_back(i);
                }
                i++;
            }
            i = nLocalInter();
            for (AtomElement &ghost_ref :inter_ghost_list) {
                if (ghost_ref.x[1] < high && ghost_ref.x[1] >= low) {
                    sendlist.push_back(i);
                }
                i++;
            }
        } else {
            low = p_domain->meas_sub_box_upper_bounding[1] - ghostlengh;
            high = p_domain->meas_sub_box_upper_bounding[1];
            i = 0;
            for (AtomElement &inter_ref :inter_list) {
                if (inter_ref.x[1] <= high && inter_ref.x[1] > low) {
                    sendlist.push_back(i);
                }
                i++;
            }
            i = nLocalInter();
            for (AtomElement &ghost_ref :inter_ghost_list) {
                if (ghost_ref.x[1] <= high && ghost_ref.x[1] > low) {
                    sendlist.push_back(i);
                }
                i++;
            }
        }
    } else {
        if (direction == 0) {
            low = p_domain->meas_sub_box_lower_bounding[2];
            high = p_domain->meas_sub_box_lower_bounding[2] + ghostlengh;
            i = 0;
            for (AtomElement &inter_ref :inter_list) {
                if (inter_ref.x[2] < high && inter_ref.x[2] >= low) {
                    sendlist.push_back(i);
                }
                i++;
            }
            i = nLocalInter();
            for (AtomElement &ghost_ref :inter_ghost_list) {
                if (ghost_ref.x[2] < high && ghost_ref.x[2] >= low) {
                    sendlist.push_back(i);
                }
                i++;
            }
        } else {
            low = p_domain->meas_sub_box_upper_bounding[2] - ghostlengh;
            high = p_domain->meas_sub_box_upper_bounding[2];
            i = 0;
            for (AtomElement &inter_ref :inter_list) {
                if (inter_ref.x[2] <= high && inter_ref.x[2] > low) {
                    sendlist.push_back(i);
                }
                i++;
            }
            i = nLocalInter();
            for (AtomElement &ghost_ref :inter_ghost_list) {
                if (ghost_ref.x[2] <= high && ghost_ref.x[2] > low) {
                    sendlist.push_back(i);
                }
                i++;
            }
        }
    }
}

unsigned long InterAtomList::getintersendnum(Domain *p_domain, int dimension, int direction) {
    interbuf.clear();
    unsigned long i = 0;
    for (AtomElement &inter_ref :inter_list) {
        if (direction == 0) {
            // we assume that, a atom cannot cross 2 or more than 2 sub-boxes
            if (inter_ref.x[dimension] < p_domain->meas_sub_box_lower_bounding[dimension]) {
                interbuf.push_back(i);
            }
        } else {
            if (inter_ref.x[dimension] >= p_domain->meas_sub_box_upper_bounding[dimension]) {
                interbuf.push_back(i);
            }
        }
        i++;
    }
    return interbuf.size();
}
