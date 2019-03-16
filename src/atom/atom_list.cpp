//
// Created by genshen on 5/8/18.
//

#include <utils/mpi_domain.h>
#include "atom_list.h"
#include "../pack/pack.h"
#include "../utils/mpi_data_types.h"

AtomList::AtomList(_type_atom_count size_x, _type_atom_count size_y, _type_atom_count size_z,
                   _type_atom_count size_sub_box_x, _type_atom_count size_sub_box_y, _type_atom_count size_sub_box_z,
                   _type_atom_count ghost_count_x, _type_atom_count ghost_count_y, _type_atom_count ghost_count_z) :
        _size(size_x * size_y * size_z),
        _size_x(size_x), _size_y(size_y), _size_z(size_z),
        _size_sub_box_x(size_sub_box_x),
        _size_sub_box_y(size_sub_box_y),
        _size_sub_box_z(size_sub_box_z),
        purge_ghost_count_x(ghost_count_x),
        purge_ghost_count_y(ghost_count_y),
        purge_ghost_count_z(ghost_count_z) {

    _atoms = new AtomElement **[size_z];
    for (_type_atom_count z = 0; z < size_z; z++) {
        _atoms[z] = new AtomElement *[size_y];
        for (_type_atom_count y = 0; y < size_y; y++) {
            _atoms[z][y] = new AtomElement[size_x];
        }
    }
}

AtomList::~AtomList() {
    for (_type_atom_count z = 0; z < _size_z; z++) {
        for (_type_atom_count y = 0; y < _size_y; y++) {
            delete[] _atoms[z][y];
        }
        delete[] _atoms[z];
    }
    delete[] _atoms;
}

void AtomList::exchangeAtomFirst(comm::Domain *p_domain, int cutlattice) {
    sendlist.resize(6);
    recvlist.resize(6);

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
    for (unsigned short d = 0; d < DIMENSION; d++) {
        for (direction = LOWER; direction <= HIGHER; direction++) {
            // 找到要发送给邻居的原子
            switch (d) {
                case 0:
                    getatomx(p_domain, cutlattice, direction, sendlist);
                    break;
                case 1:
                    getatomy(p_domain, cutlattice, direction, sendlist);
                    break;
                case 2:
                    getatomz(p_domain, cutlattice, direction, sendlist);
                    break;
                default:
                    break;
            }

            // 当原子要跨越周期性边界, 原子坐标必须要做出调整
            // only one periodic boundary appears at one dimension, so the length of array can be 3, not 6.
            double offset[DIMENSION] = {0.0, 0.0, 0.0};
            // 进程在左侧边界
            if (p_domain->grid_coord_sub_box[d] == 0 && direction == LOWER) {
                offset[d] = p_domain->meas_global_length[d];
            }

            // 进程在右侧边界
            if (p_domain->grid_coord_sub_box[d] == p_domain->grid_size[d] - 1 && direction == HIGHER) {
                offset[d] = -((p_domain->meas_global_length[d]));
            }

            // 初始化发送缓冲区
            numPartsToSend[d][direction] = sendlist[iswap].size();
            sendbuf[direction] = new LatParticleData[numPartsToSend[d][direction]];
            pack::pack_send(d, numPartsToSend[d][direction], offset, *this,
                            sendbuf[direction], sendlist[iswap++]);
//            _atom->pack_send(d, numPartsToSend[d][direction], sendlist[iswap++], sendbuf[direction], offset);
        }

        // 与上下邻居通信
        // send ghost atoms.
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numsend = numPartsToSend[d][direction];
            int numrecv;

            // 向下/上发送并从上/下接收
            MPI_Isend(sendbuf[direction], numsend, mpi_types::_mpi_latParticle_data,
                      p_domain->rank_id_neighbours[d][direction],
                      99,
                      MPIDomain::sim_processor.comm, &send_requests[d][direction]);
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
        // receive ghost atoms
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numrecv = numPartsToRecv[d][direction];
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);

            //将收到的粒子位置信息加到对应存储位置上
            pack::unpack_recvfirst(d, direction, numrecv, *this,
                                   p_domain->dbx_lattice_size_ghost,
                                   p_domain->dbx_lattice_size_sub_box,
                                   p_domain->dbx_lattice_size_ghost_extended,
                                   recvbuf[direction], recvlist);
//            _atom->unpack_recvfirst(d, direction, numrecv, recvbuf[direction], recvlist);

            // 释放buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}

void AtomList::exchangeAtom(comm::Domain *p_domain) {
    double ghostlengh[DIMENSION]; // ghost区域大小

    for (int d = 0; d < DIMENSION; d++) {
        ghostlengh[d] = p_domain->meas_ghost_length[d];
    }

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
    for (unsigned short d = 0; d < DIMENSION; d++) {
        for (direction = LOWER; direction <= HIGHER; direction++) {
            // 当原子要跨越周期性边界, 原子坐标必须要做出调整
            double offset[DIMENSION] = {0.0, 0.0, 0.0};
            // 进程在左侧边界
            if (p_domain->grid_coord_sub_box[d] == 0 && direction == LOWER) {
                offset[d] = p_domain->meas_global_length[d];
            }
            // 进程在右侧边界
            if (p_domain->grid_coord_sub_box[d] == p_domain->grid_size[d] - 1 && direction == HIGHER) {
                offset[d] = -((p_domain->meas_global_length[d]));
            }

            // 初始化发送缓冲区
            numPartsToSend[d][direction] = sendlist[iswap].size();
            sendbuf[direction] = new LatParticleData[numPartsToSend[d][direction]];
            pack::pack_send(d, numPartsToSend[d][direction], offset, *this,
                            sendbuf[direction], sendlist[iswap++]);
//            _atom->pack_send(d, numPartsToSend[d][direction], sendlist[iswap++], sendbuf[direction], offset);
        }

        // 与上下邻居通信
        for (direction = LOWER; direction <= HIGHER; direction++) {

            int numsend = numPartsToSend[d][direction];
            int numrecv;

            MPI_Isend(sendbuf[direction], numsend, mpi_types::_mpi_latParticle_data,
                      p_domain->rank_id_neighbours[d][direction],
                      99, MPIDomain::sim_processor.comm,
                      &send_requests[d][direction]);
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
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);

            //将收到的粒子位置信息加到对应存储位置上
            pack::unpack_recv(d, direction, numrecv, *this, recvbuf[direction], recvlist);
//            _atom->unpack_recv(d, direction, numrecv, recvbuf[direction], recvlist);

            // 释放buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}


void AtomList::getatomx(comm::Domain *p_domain, int _cutlattice, int direction,
                        std::vector<std::vector<_type_atom_id> > &sendlist) {
    _type_atom_id i;
    if (direction == 0) {
        //找到要发送到邻居进程的区域
        int xstart = p_domain->dbx_lattice_size_ghost[0];
        int ystart = p_domain->dbx_lattice_size_ghost[1];
        int zstart = p_domain->dbx_lattice_size_ghost[2];
        int xstop = xstart + p_domain->dbx_lattice_size_ghost[0]; // note: this is ghost lattice size.
        int ystop = ystart + p_domain->dbx_lattice_size_sub_box[1];
        int zstop = zstart + p_domain->dbx_lattice_size_sub_box[2];

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[0].push_back(i);
                }
            }
        }
    } else {
        //找到要发送到邻居进程的区域
        int xstart = p_domain->dbx_lattice_size_ghost[0] + p_domain->dbx_lattice_size_sub_box[0] - ((_cutlattice) * 2);
        int ystart = p_domain->dbx_lattice_size_ghost[1];
        int zstart = p_domain->dbx_lattice_size_ghost[2];
        int xstop = p_domain->dbx_lattice_size_ghost[0] + p_domain->dbx_lattice_size_sub_box[0];
        int ystop = ystart + p_domain->dbx_lattice_size_sub_box[1];
        int zstop = zstart + p_domain->dbx_lattice_size_sub_box[2];

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[1].push_back(i);
                }
            }
        }
    }
}

void AtomList::getatomy(comm::Domain *p_domain, int _cutlattice, int direction,
                        std::vector<std::vector<_type_atom_id> > &sendlist) {
    int i;
    if (direction == 0) {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = p_domain->dbx_lattice_size_ghost[1];
        int zstart = p_domain->dbx_lattice_size_ghost[2];
        int xstop = p_domain->dbx_lattice_size_ghost_extended[0];
        int ystop = ystart + _cutlattice;
        int zstop = zstart + p_domain->dbx_lattice_size_sub_box[2];

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[2].push_back(i);
                }
            }
        }
    } else {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = p_domain->dbx_lattice_size_ghost[1] + p_domain->dbx_lattice_size_sub_box[1] - (_cutlattice);
        int zstart = p_domain->dbx_lattice_size_ghost[2];
        int xstop = p_domain->dbx_lattice_size_ghost_extended[0];
        int ystop = p_domain->dbx_lattice_size_ghost[1] + p_domain->dbx_lattice_size_sub_box[1];
        int zstop = zstart + p_domain->dbx_lattice_size_sub_box[2];

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[3].push_back(i);
                }
            }
        }
    }
}

void AtomList::getatomz(comm::Domain *p_domain, int _cutlattice, int direction,
                        std::vector<std::vector<_type_atom_id> > &sendlist) {
    int i;
    if (direction == 0) {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = 0;
        int zstart = p_domain->dbx_lattice_size_ghost[2];
        int xstop = p_domain->dbx_lattice_size_ghost_extended[0];
        int ystop = p_domain->dbx_lattice_size_ghost_extended[1];
        int zstop = zstart + _cutlattice;

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[4].push_back(i);
                }
            }
        }
    } else {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = 0;
        int zstart = p_domain->dbx_lattice_size_ghost[2] + p_domain->dbx_lattice_size_sub_box[2] - (_cutlattice);
        int xstop = p_domain->dbx_lattice_size_ghost_extended[0];
        int ystop = p_domain->dbx_lattice_size_ghost_extended[1];
        int zstop = p_domain->dbx_lattice_size_ghost[2] + p_domain->dbx_lattice_size_sub_box[2];

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[5].push_back(i);
                }
            }
        }
    }
}

/**
 *******************************************
 * following (to end of this file) are the
 * implementations about atom list iterator.
 *******************************************
 */
AtomList::iterator AtomList::begin() {
    return AtomList::iterator(nullptr);
}

AtomList::iterator AtomList::end() {
    return AtomList::iterator(nullptr);
}

AtomList::const_iterator AtomList::begin() const {
    return AtomList::const_iterator(nullptr);
}

AtomList::const_iterator AtomList::end() const {
    return AtomList::const_iterator(nullptr);
}

// iterator
AtomList::iterator::iterator(AtomElement *ptr) : ptr_(ptr) {}

AtomList::iterator::self_type &AtomList::iterator::operator++() {
    ptr_++;
    return *this;
}

AtomList::iterator::self_type AtomList::iterator::operator++(int) {
    self_type i = *this;
    ptr_++;
    return i;
}

AtomElement &AtomList::iterator::operator*() {
    return *ptr_;
}

AtomElement *AtomList::iterator::operator->() {
    return ptr_;
}

bool AtomList::iterator::operator==(const AtomList::iterator::self_type &rhs) {
    return ptr_ == rhs.ptr_;
}

bool AtomList::iterator::operator!=(const AtomList::iterator::self_type &rhs) {
    return ptr_ != rhs.ptr_;
}

void AtomList::iterator::next() {

}

// const iterator
AtomList::const_iterator::const_iterator(AtomElement *ptr) : ptr_(ptr) {}

AtomList::const_iterator::self_type &AtomList::const_iterator::operator++() {
    ptr_++;
    return *this;
}

AtomList::const_iterator::self_type AtomList::const_iterator::operator++(int junk) {
    self_type i = *this;
    ptr_++;
    return i;
}

const AtomElement &AtomList::const_iterator::operator*() {
    return *ptr_;
}

bool AtomList::const_iterator::operator==(const AtomList::const_iterator::self_type &rhs) {
    return ptr_ == rhs.ptr_;
}

bool AtomList::const_iterator::operator!=(const AtomList::const_iterator::self_type &rhs) {
    return ptr_ != rhs.ptr_;
}
