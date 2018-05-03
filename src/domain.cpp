#include <iostream>
#include <logs/logs.h>
#include "domain.h"

DomainDecomposition::DomainDecomposition() {
    // 3维拓扑
    for (int d = 0; d < DIMENSION; d++) {
        _grid_size[d] = 0;
        _globalLength[d] = 0;
    }
}

DomainDecomposition::~DomainDecomposition() {
    MPI_Type_free(&_mpi_Particle_data);
    MPI_Type_free(&_mpi_latParticle_data);
}

DomainDecomposition *DomainDecomposition::decomposition() {
    // Assume N can be decomposed as N = N_x * N_y * N_z,
    // then we have: _grid_size[0] = N_x, _grid_size[1] = N_y, _grid_size[1] = N_z.
    // Fill in the _grid_size array such that the product of _grid_size[i] for i=0 to DIMENSION-1 equals N.
    MPI_Dims_create(kiwi::mpiUtils::all_ranks, DIMENSION, (int *) &_grid_size);
    if (kiwi::mpiUtils::own_rank == MASTER_PROCESSOR) {
        kiwi::logs::i("decomposition", "MPI grid dimensions: {0},{1},{2}\n",
                      _grid_size[0], _grid_size[1], _grid_size[2]);
    }

    int period[DIMENSION];
    // 3维拓扑
    for (int d = 0; d < DIMENSION; d++) {
        period[d] = 1;
    }
    // sort the processors to fit 3D cartesian topology.
    // the rank id may change.
    MPI_Comm _comm;
    int _debug_old_rank = kiwi::mpiUtils::own_rank;
    MPI_Cart_create(MPI_COMM_WORLD, DIMENSION, _grid_size, period, true, &_comm);
    kiwi::mpiUtils::onGlobalCommChanged(_comm);

    // get cartesian coordinate of current processor.
    MPI_Cart_coords(kiwi::mpiUtils::global_comm, kiwi::mpiUtils::own_rank, DIMENSION, _coords);
    kiwi::logs::d("decomposition", "old_rank_id: {0}, MPI coordinate of current process: x:{1},y{2},z{3}\n",
                  _debug_old_rank, _coords[0], _coords[1], _coords[2]);

    // get the rank ids of contiguous processors of current processor.
    for (int d = 0; d < DIMENSION; d++) {
        MPI_Cart_shift(_comm, d, 1, &_rank_id_neighbours[d][LOWER], &_rank_id_neighbours[d][HIGHER]);
    }

    particledata::setMPIType(_mpi_Particle_data);  // todo move code to other place?
    LatParticleData::setMPIType(_mpi_latParticle_data);
    intersendlist.resize(6);
    interrecvlist.resize(6);
    return this;
}

DomainDecomposition *DomainDecomposition::createGlobalDomain(const int64_t *phaseSpace, const double latticeConst) {
    for (int d = 0; d < DIMENSION; d++) {
        //phaseSpace个单位长度(单位长度即latticeconst)
        _globalLength[d] = phaseSpace[d] * latticeConst;
        _coord_global_box_low[d] = 0;
        _coord_global_box_high[d] = _globalLength[d];
    }
    return this;
}

DomainDecomposition *DomainDecomposition::createLocalBoxDomain(const int64_t *phaseSpace,
                                                               const double latticeConst, const double cutoffRadius) {
    for (int d = 0; d < DIMENSION; d++) {
        _boundingBoxMin[d] = getBoundingBoxMin(d);
        _boundingBoxMax[d] = getBoundingBoxMax(d);

        _ghostLength[d] = cutoffRadius; // todo ??

        _ghostBoundingBoxMin[d] = _boundingBoxMin[d] - _ghostLength[d];
        _ghostBoundingBoxMax[d] = _boundingBoxMax[d] + _ghostLength[d];
    }
    return this;
}

double DomainDecomposition::getSubBoxLowerBounding(int dimension) const { // todo inline.
    return _boundingBoxMin[dimension];
}

double DomainDecomposition::getUpperBoxLowerBounding(int dimension) const {
    return _boundingBoxMax[dimension];
}

double DomainDecomposition::getGhostLength(int index) const {
    return _ghostLength[index];
}

double DomainDecomposition::getBoundingBoxMin(int dimension) {
    return _coords[dimension] * (getGlobalLength(dimension) / _grid_size[dimension]);
}

double DomainDecomposition::getBoundingBoxMax(int dimension) {
    return (_coords[dimension] + 1) * (getGlobalLength(dimension) / _grid_size[dimension]);
}

void DomainDecomposition::exchangeAtomfirst(atom *_atom) {

    double ghostlengh[DIMENSION]; // ghost区域大小

    for (int d = 0; d < DIMENSION; d++) {
        ghostlengh[d] = getGhostLength(d);
    }
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
        // 当原子要跨越周期性边界, 原子坐标必须要做出调整

        double offsetLower[DIMENSION];
        double offsetHigher[DIMENSION];
        offsetLower[d] = 0.0;
        offsetHigher[d] = 0.0;

        // 进程在左侧边界
        if (_coords[d] == 0)
            offsetLower[d] = getGlobalLength(d);
        // 进程在右侧边界
        if (_coords[d] == _grid_size[d] - 1)
            offsetHigher[d] = -(getGlobalLength(d));

        for (direction = LOWER; direction <= HIGHER; direction++) {
            // 找到要发送给邻居的原子
            if (d == 0) {
                _atom->getatomx(direction, sendlist);
            } else if (d == 1) {
                _atom->getatomy(direction, sendlist);
            } else {
                _atom->getatomz(direction, sendlist);
            }

            double shift = 0.0;
            if (direction == LOWER)
                shift = offsetLower[d];
            if (direction == HIGHER)
                shift = offsetHigher[d];

            // 初始化发送缓冲区
            numPartsToSend[d][direction] = sendlist[iswap].size();
            sendbuf[direction] = new LatParticleData[numPartsToSend[d][direction]];
            _atom->pack_send(d, numPartsToSend[d][direction], sendlist[iswap++], sendbuf[direction], shift);
        }

        // 与上下邻居通信
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numsend = numPartsToSend[d][direction];
            int numrecv;

            // 向下/上发送并从上/下接收
            MPI_Isend(sendbuf[direction], numsend, _mpi_latParticle_data, _rank_id_neighbours[d][direction], 99,
                      kiwi::mpiUtils::global_comm,
                      &send_requests[d][direction]);
            MPI_Probe(_rank_id_neighbours[d][(direction + 1) % 2], 99, kiwi::mpiUtils::global_comm,
                      &status);//测试邻居是否有信息发送给本地
            MPI_Get_count(&status, _mpi_latParticle_data, &numrecv);//得到要接收的粒子数目
            // 初始化接收缓冲区
            //依据得到发送方要发送粒子信息大小，初始化接收缓冲区
            recvbuf[direction] = new LatParticleData[numrecv];
            numPartsToRecv[d][direction] = numrecv;
            MPI_Irecv(recvbuf[direction], numrecv, _mpi_latParticle_data, _rank_id_neighbours[d][(direction + 1) % 2],
                      99,
                      kiwi::mpiUtils::global_comm, &recv_requests[d][direction]);
        }

        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numrecv = numPartsToRecv[d][direction];
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);

            //将收到的粒子位置信息加到对应存储位置上
            _atom->unpack_recvfirst(d, direction, numrecv, recvbuf[direction], recvlist);

            // 释放buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}

void DomainDecomposition::exchangeAtom(atom *_atom) {

    double ghostlengh[DIMENSION]; // ghost区域大小

    for (int d = 0; d < DIMENSION; d++) {
        ghostlengh[d] = getGhostLength(d);
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
        // 当原子要跨越周期性边界, 原子坐标必须要做出调整

        double offsetLower[DIMENSION];
        double offsetHigher[DIMENSION];
        offsetLower[d] = 0.0;
        offsetHigher[d] = 0.0;

        // 进程在左侧边界
        if (_coords[d] == 0)
            offsetLower[d] = getGlobalLength(d);
        // 进程在右侧边界
        if (_coords[d] == _grid_size[d] - 1)
            offsetHigher[d] = -(getGlobalLength(d));

        for (direction = LOWER; direction <= HIGHER; direction++) {
            double shift = 0.0;
            if (direction == LOWER)
                shift = offsetLower[d];
            if (direction == HIGHER)
                shift = offsetHigher[d];

            // 初始化发送缓冲区
            numPartsToSend[d][direction] = sendlist[iswap].size();
            sendbuf[direction] = new LatParticleData[numPartsToSend[d][direction]];
            _atom->pack_send(d, numPartsToSend[d][direction], sendlist[iswap++], sendbuf[direction], shift);
        }

        // 与上下邻居通信
        for (direction = LOWER; direction <= HIGHER; direction++) {

            int numsend = numPartsToSend[d][direction];
            int numrecv;

            MPI_Isend(sendbuf[direction], numsend, _mpi_latParticle_data, _rank_id_neighbours[d][direction], 99,
                      kiwi::mpiUtils::global_comm,
                      &send_requests[d][direction]);
            MPI_Probe(_rank_id_neighbours[d][(direction + 1) % 2], 99, kiwi::mpiUtils::global_comm,
                      &status);//测试邻居是否有信息发送给本地
            MPI_Get_count(&status, _mpi_latParticle_data, &numrecv);//得到要接收的粒子数目
            // 初始化接收缓冲区
            //依据得到发送方要发送粒子信息大小，初始化接收缓冲区
            recvbuf[direction] = new LatParticleData[numrecv];
            numPartsToRecv[d][direction] = numrecv;
            MPI_Irecv(recvbuf[direction], numrecv, _mpi_latParticle_data, _rank_id_neighbours[d][(direction + 1) % 2],
                      99,
                      kiwi::mpiUtils::global_comm, &recv_requests[d][direction]);
        }

        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numrecv = numPartsToRecv[d][direction];
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);

            //将收到的粒子位置信息加到对应存储位置上
            _atom->unpack_recv(d, direction, numrecv, recvbuf[direction], recvlist);

            // 释放buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}

void DomainDecomposition::exchangeInter(atom *_atom) {
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
            numPartsToSend[d][direction] = _atom->getintersendnum(d, direction);
            sendbuf[direction] = new particledata[numPartsToSend[d][direction]];

            _atom->pack_intersend(sendbuf[direction]);
        }

        // 与上下邻居通信
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numsend = numPartsToSend[d][direction];
            int numrecv;

            MPI_Isend(sendbuf[direction], numsend, _mpi_Particle_data, _rank_id_neighbours[d][direction], 99,
                      kiwi::mpiUtils::global_comm, &send_requests[d][direction]);
            MPI_Probe(_rank_id_neighbours[d][(direction + 1) % 2], 99, kiwi::mpiUtils::global_comm,
                      &status);//测试邻居是否有信息发送给本地
            MPI_Get_count(&status, _mpi_Particle_data, &numrecv);//得到要接收的粒子数目
            // 初始化接收缓冲区
            //依据得到发送方要发送粒子信息大小，初始化接收缓冲区
            recvbuf[direction] = new particledata[numrecv];
            numPartsToRecv[d][direction] = numrecv;
            MPI_Irecv(recvbuf[direction], numrecv, _mpi_Particle_data, _rank_id_neighbours[d][(direction + 1) % 2], 99,
                      kiwi::mpiUtils::global_comm, &recv_requests[d][direction]);
        }

        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numrecv = numPartsToRecv[d][direction];
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);

            //将收到的粒子位置信息加到对应存储位置上
            _atom->unpack_interrecv(d, numrecv, recvbuf[direction]);

            // 释放buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}

void DomainDecomposition::borderInter(atom *_atom) {
    double ghostlengh[DIMENSION]; // ghost区域大小

    for (int d = 0; d < DIMENSION; d++) {
        ghostlengh[d] = getGhostLength(d);
    }

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
        if (_coords[d] == 0)
            offsetLower[d] = getGlobalLength(d);
        // 进程在右侧边界
        if (_coords[d] == _grid_size[d] - 1)
            offsetHigher[d] = -(getGlobalLength(d));

        for (direction = LOWER; direction <= HIGHER; direction++) {
            // 找到要发送给邻居的原子
            _atom->getIntertosend(d, direction, ghostlengh[d], intersendlist[iswap]);

            double shift = 0.0;
            if (direction == LOWER)
                shift = offsetLower[d];
            if (direction == HIGHER)
                shift = offsetHigher[d];

            // 初始化发送缓冲区
            numPartsToSend[d][direction] = intersendlist[iswap].size();
            sendbuf[direction] = new LatParticleData[numPartsToSend[d][direction]];
            _atom->pack_bordersend(d, numPartsToSend[d][direction], intersendlist[iswap++], sendbuf[direction], shift);

        }

        // 与上下邻居通信
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numsend = numPartsToSend[d][direction];
            int numrecv;

            MPI_Isend(sendbuf[direction], numsend, _mpi_latParticle_data, _rank_id_neighbours[d][direction], 99,
                      kiwi::mpiUtils::global_comm, &send_requests[d][direction]);
            MPI_Probe(_rank_id_neighbours[d][(direction + 1) % 2], 99, kiwi::mpiUtils::global_comm,
                      &status);//测试邻居是否有信息发送给本地
            MPI_Get_count(&status, _mpi_latParticle_data, &numrecv);//得到要接收的粒子数目
            // 初始化接收缓冲区
            //依据得到发送方要发送粒子信息大小，初始化接收缓冲区
            recvbuf[direction] = new LatParticleData[numrecv];
            numPartsToRecv[d][direction] = numrecv;
            MPI_Irecv(recvbuf[direction], numrecv, _mpi_latParticle_data, _rank_id_neighbours[d][(direction + 1) % 2],
                      99,
                      kiwi::mpiUtils::global_comm, &recv_requests[d][direction]);
        }

        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numrecv = numPartsToRecv[d][direction];
            interrecvlist[jswap].resize(numrecv);
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);

            //将收到的粒子位置信息加到对应存储位置上
            _atom->unpack_borderrecv(numrecv, recvbuf[direction], interrecvlist[jswap++]);

            // 释放buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}

void DomainDecomposition::sendrho(atom *_atom) {
    // 发送、接收数据缓冲区
    int numPartsToSend[DIMENSION][2];
    int numPartsToRecv[DIMENSION][2];
    double *sendbuf[2];
    double *recvbuf[2];

    MPI_Status status;
    MPI_Status send_statuses[DIMENSION][2];
    MPI_Status recv_statuses[DIMENSION][2];
    MPI_Request send_requests[DIMENSION][2];
    MPI_Request recv_requests[DIMENSION][2];

    int direction;
    int iswap = 5;
    for (int d = (DIMENSION - 1); d >= 0; d--) {
        for (direction = LOWER; direction <= HIGHER; direction++) {
            numPartsToSend[d][direction] = recvlist[iswap].size();
            sendbuf[direction] = new double[numPartsToSend[d][direction]];
            _atom->pack_rho(numPartsToSend[d][direction], recvlist[iswap--], sendbuf[direction]);
        }
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numsend = numPartsToSend[d][direction];
            int numrecv;

            MPI_Isend(sendbuf[direction], numsend, MPI_DOUBLE, _rank_id_neighbours[d][direction], 99,
                      kiwi::mpiUtils::global_comm, &send_requests[d][direction]);
            MPI_Probe(_rank_id_neighbours[d][(direction + 1) % 2], 99, kiwi::mpiUtils::global_comm,
                      &status);//测试邻居是否有信息发送给本地
            MPI_Get_count(&status, MPI_DOUBLE, &numrecv);//得到要接收的粒子数目
            // 初始化接收缓冲区
            //依据得到发送方要发送粒子信息大小，初始化接收缓冲区
            recvbuf[direction] = new double[numrecv];
            numPartsToRecv[d][direction] = numrecv;
            MPI_Irecv(recvbuf[direction], numrecv, MPI_DOUBLE, _rank_id_neighbours[d][(direction + 1) % 2], 99,
                      kiwi::mpiUtils::global_comm, &recv_requests[d][direction]);
        }
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numrecv = numPartsToRecv[d][direction];
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);

            //将收到的电子云密度信息加到对应存储位置上
            _atom->unpack_rho(d, direction, recvbuf[direction], sendlist);

            // 释放buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}

void DomainDecomposition::sendDfEmbed(atom *_atom) {
    // 发送、接收数据缓冲区
    int numPartsToSend[DIMENSION][2];
    int numPartsToRecv[DIMENSION][2];
    double *sendbuf[2];
    double *recvbuf[2];

    MPI_Status status;
    MPI_Status send_statuses[DIMENSION][2];
    MPI_Status recv_statuses[DIMENSION][2];
    MPI_Request send_requests[DIMENSION][2];
    MPI_Request recv_requests[DIMENSION][2];

    int direction;
    int iswap = 0;
    int jswap = 0;
    for (unsigned short d = 0; d < DIMENSION; d++) {
        for (direction = LOWER; direction <= HIGHER; direction++) {
            // 初始化发送缓冲区
            numPartsToSend[d][direction] = sendlist[iswap].size() + intersendlist[iswap].size();
            sendbuf[direction] = new double[numPartsToSend[d][direction]];
            _atom->pack_df(sendlist[iswap], intersendlist[iswap], sendbuf[direction]);
            iswap++;
        }

        // 与上下邻居通信
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numsend = numPartsToSend[d][direction];
            int numrecv;

            MPI_Isend(sendbuf[direction], numsend, MPI_DOUBLE, _rank_id_neighbours[d][direction], 99,
                      kiwi::mpiUtils::global_comm, &send_requests[d][direction]);
            MPI_Probe(_rank_id_neighbours[d][(direction + 1) % 2], 99, kiwi::mpiUtils::global_comm,
                      &status);//测试邻居是否有信息发送给本地
            MPI_Get_count(&status, MPI_DOUBLE, &numrecv);//得到要接收的粒子数目
            // 初始化接收缓冲区
            //依据得到发送方要发送粒子信息大小，初始化接收缓冲区
            recvbuf[direction] = new double[numrecv];
            numPartsToRecv[d][direction] = numrecv;
            MPI_Irecv(recvbuf[direction], numrecv, MPI_DOUBLE, _rank_id_neighbours[d][(direction + 1) % 2], 99,
                      kiwi::mpiUtils::global_comm, &recv_requests[d][direction]);
        }

        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numrecv = numPartsToRecv[d][direction];
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);

            //将收到的嵌入能导数信息加到对应存储位置上
            _atom->unpack_df(numrecv, recvbuf[direction], recvlist[jswap], interrecvlist[jswap]);
            jswap++;

            // 释放buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}

void DomainDecomposition::sendForce(atom *_atom) {
    // 发送、接收数据缓冲区
    int numPartsToSend[DIMENSION][2];
    int numPartsToRecv[DIMENSION][2];
    double *sendbuf[2];
    double *recvbuf[2];

    MPI_Status status;
    MPI_Status send_statuses[DIMENSION][2];
    MPI_Status recv_statuses[DIMENSION][2];
    MPI_Request send_requests[DIMENSION][2];
    MPI_Request recv_requests[DIMENSION][2];

    int direction;
    int iswap = 5;
    for (int d = (DIMENSION - 1); d >= 0; d--) {
        for (direction = LOWER; direction <= HIGHER; direction++) {
            numPartsToSend[d][direction] = recvlist[iswap].size();
            sendbuf[direction] = new double[numPartsToSend[d][direction] * 3];
            _atom->pack_force(numPartsToSend[d][direction], recvlist[iswap--], sendbuf[direction]);
        }
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numsend = numPartsToSend[d][direction] * 3;
            int numrecv;

            MPI_Isend(sendbuf[direction], numsend, MPI_DOUBLE, _rank_id_neighbours[d][direction], 99,
                      kiwi::mpiUtils::global_comm, &send_requests[d][direction]);
            MPI_Probe(_rank_id_neighbours[d][(direction + 1) % 2], 99, kiwi::mpiUtils::global_comm,
                      &status);//测试邻居是否有信息发送给本地
            MPI_Get_count(&status, MPI_DOUBLE, &numrecv);//得到要接收的粒子数目
            // 初始化接收缓冲区
            //依据得到发送方要发送粒子信息大小，初始化接收缓冲区
            recvbuf[direction] = new double[numrecv];
            numPartsToRecv[d][direction] = numrecv;
            MPI_Irecv(recvbuf[direction], numrecv, MPI_DOUBLE, _rank_id_neighbours[d][(direction + 1) % 2], 99,
                      kiwi::mpiUtils::global_comm, &recv_requests[d][direction]);
        }
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numrecv = numPartsToRecv[d][direction]; // todo remove not used variable.
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);

            //将收到的粒子位置信息加到对应存储位置上
            _atom->unpack_force(d, direction, recvbuf[direction], sendlist);

            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}
