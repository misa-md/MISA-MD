#include <logs/logs.h>
#include "domain.h"
#include "utils/mpi_domain.h"
#include "utils/mpi_data_types.h"
#include "pack/pack.h"

Domain::Domain(const int64_t *phaseSpace, const double latticeConst,
               const double cutoffRadiusFactor) :
        _lattice_const(latticeConst), _cutoff_radius_factor(cutoffRadiusFactor) {
    // 3维拓扑
    for (int d = 0; d < DIMENSION; d++) {
        _phase_space[d] = phaseSpace[d]; // initialize phase space.
        _grid_size[d] = 0;
        _meas_global_length[d] = 0;
    }
}

Domain::~Domain() {}

Domain *Domain::decomposition() {
    // Assume N can be decomposed as N = N_x * N_y * N_z,
    // then we have: _grid_size[0] = N_x, _grid_size[1] = N_y, _grid_size[1] = N_z.
    // Fill in the _grid_size array such that the product of _grid_size[i] for i=0 to DIMENSION-1 equals N.
    MPI_Dims_create(MPIDomain::sim_processor.all_ranks, DIMENSION,
                    _grid_size); // fixme origin code: (int *) &_grid_size
    kiwi::logs::i(MASTER_PROCESSOR, "decomposition", "MPI grid dimensions: {0},{1},{2}\n",
                  _grid_size[0], _grid_size[1], _grid_size[2]);

    int period[DIMENSION];
    // 3维拓扑
    for (int d = 0; d < DIMENSION; d++) {
        period[d] = 1;
    }
    // sort the processors to fit 3D cartesian topology.
    // the rank id may change.
    MPI_Comm _comm;
    int _debug_old_rank = MPIDomain::sim_processor.own_rank;
    MPI_Cart_create(MPI_COMM_WORLD, DIMENSION, _grid_size, period, true, &_comm);
    kiwi::mpiUtils::onGlobalCommChanged(_comm);
    MPIDomain::sim_processor = kiwi::mpiUtils::global_process; // set new domain.

    // get cartesian coordinate of current processor.
    MPI_Cart_coords(MPIDomain::sim_processor.comm, MPIDomain::sim_processor.own_rank, DIMENSION, _grid_coord_sub_box);
    kiwi::logs::d("decomposition", "old_rank_id: {0}, MPI coordinate of current process: x:{1},y{2},z{3}\n",
                  _debug_old_rank, _grid_coord_sub_box[0], _grid_coord_sub_box[1], _grid_coord_sub_box[2]);

    // get the rank ids of contiguous processors of current processor.
    for (int d = 0; d < DIMENSION; d++) {
        MPI_Cart_shift(_comm, d, 1, &_rank_id_neighbours[d][LOWER], &_rank_id_neighbours[d][HIGHER]);
    }
    return this;
}

Domain *Domain::createGlobalDomain() {
    for (int d = 0; d < DIMENSION; d++) {
        //phaseSpace个单位长度(单位长度即latticeconst)
        _meas_global_length[d] = _phase_space[d] * _lattice_const;
        _meas_global_box_coord_lower[d] = 0; // lower bounding is set to 0 by default.
        _meas_global_box_coord_upper[d] = _meas_global_length[d];
    }
    return this;
}

Domain *Domain::createSubBoxDomain() {
    for (int d = 0; d < DIMENSION; d++) {
        // the lower and upper bounding of current sub-box.
        _meas_sub_box_lower_bounding[d] = _meas_global_box_coord_lower[d] +
                                          _grid_coord_sub_box[d] * (_meas_global_length[d] / _grid_size[d]);
        _meas_sub_box_upper_bounding[d] = _meas_global_box_coord_lower[d] +
                                          (_grid_coord_sub_box[d] + 1) * (_meas_global_length[d] / _grid_size[d]);

        _meas_ghost_length[d] = _cutoff_radius_factor * _lattice_const; // ghost length todo

        _meas_ghost_lower_bounding[d] = _meas_sub_box_lower_bounding[d] - _meas_ghost_length[d];
        _meas_ghost_upper_bounding[d] = _meas_sub_box_upper_bounding[d] + _meas_ghost_length[d];
    }

    // set lattice size of local sub-box.
    for (int d = 0; d < DIMENSION; d++) {
        _lattice_size_sub_box[d] = (_grid_coord_sub_box[d] + 1) * _phase_space[d] / _grid_size[d] -
                                   (_grid_coord_sub_box[d]) * _phase_space[d] / _grid_size[d];
    }
    _lattice_size_sub_box[0] *= 2; // todo ?? why

    /*
    nghostx = p_domain->getSubBoxLatticeSize(0) + 2 * 2 * ( ceil( cutoffRadius / _latticeconst ) + 1 );
    nghosty = p_domain->getSubBoxLatticeSize(1) + 2 * ( ceil( cutoffRadius / _latticeconst ) + 1 );
    nghostz = p_domain->getSubBoxLatticeSize(2) + 2 * ( ceil( cutoffRadius / _latticeconst ) + 1 );
    */
    // set ghost lattice size.
    for (int d = 0; d < DIMENSION; d++) {
        // i * ceil(x) >= ceil(i*x) for all x ∈ R and i ∈ Z
        _lattice_size_ghost[d] = (d == 0) ? 2 * ceil(_cutoff_radius_factor) : ceil(_cutoff_radius_factor);
        _lattice_size_ghost_extended[d] = _lattice_size_sub_box[d] + 2 * _lattice_size_ghost[d];
    }

    // set lattice coordinate boundary in global and local coordinate system(GCS and LCS).
    setSubBoxDomainGCS();
    setSubBoxDomainLCS();
    return this;
}

void Domain::setSubBoxDomainGCS() {
    // set lattice coordinate boundary of sub-box.
    for (int d = 0; d < DIMENSION; d++) {
        // floor equals to "/" if all operation number >=0.
        _lattice_coord_sub_box_lower[d] = _grid_coord_sub_box[d] * _phase_space[d] /
                                          _grid_size[d]; // todo set measure coord = lower*lattice_const.
        _lattice_coord_sub_box_upper[d] = (_grid_coord_sub_box[d] + 1) * _phase_space[d] / _grid_size[d];
    }
    _lattice_coord_sub_box_lower[0] *= 2;
    _lattice_coord_sub_box_upper[0] *= 2;

    /*
   loghostx = p_domain->getSubBoxLatticeCoordLower(0) - 2 * ( ceil( cutoffRadius / _latticeconst ) + 1 );
   loghosty = p_domain->getSubBoxLatticeCoordLower(1) - ( ceil( cutoffRadius / _latticeconst ) + 1 );
   loghostz = p_domain->getGlobalSubBoxLatticeCoordLower(2) - ( ceil( cutoffRadius / _latticeconst ) + 1 );
   */
    // set lattice coordinate boundary for ghost.
    for (int d = 0; d < DIMENSION; d++) {
        // todo too integer minus, cut too many??
        _lattice_coord_ghost_lower[d] = _lattice_coord_sub_box_lower[d] - _lattice_size_ghost[d];
        _lattice_coord_ghost_upper[d] = _lattice_coord_sub_box_upper[d] + _lattice_size_ghost[d];
    }
}

void Domain::setSubBoxDomainLCS() {
    for (int d = 0; d < DIMENSION; d++) {
        _local_lattice_coord_ghost_lower[d] = 0;
        _local_lattice_coord_ghost_upper[d] = _lattice_size_ghost_extended[d];
        _local_lattice_coord_sub_box_lower[d] = _lattice_size_ghost[d];
        _local_lattice_coord_sub_box_upper[d] = _lattice_size_ghost[d] + _lattice_size_sub_box[d];
    }
}

void Domain::sendrho(atom *_atom) {
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
            pack::pack_rho(numPartsToSend[d][direction], _atom->getAtomListRef(),
                           sendbuf[direction], recvlist[iswap--]);
//            _atom->pack_rho(numPartsToSend[d][direction], recvlist[iswap--], sendbuf[direction]);
        }
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numsend = numPartsToSend[d][direction];
            int numrecv;

            MPI_Isend(sendbuf[direction], numsend, MPI_DOUBLE, _rank_id_neighbours[d][direction], 99,
                      MPIDomain::sim_processor.comm, &send_requests[d][direction]);
            MPI_Probe(_rank_id_neighbours[d][(direction + 1) % 2], 99, MPIDomain::sim_processor.comm,
                      &status);//测试邻居是否有信息发送给本地
            MPI_Get_count(&status, MPI_DOUBLE, &numrecv);//得到要接收的粒子数目
            // 初始化接收缓冲区
            //依据得到发送方要发送粒子信息大小，初始化接收缓冲区
            recvbuf[direction] = new double[numrecv];
            numPartsToRecv[d][direction] = numrecv;
            MPI_Irecv(recvbuf[direction], numrecv, MPI_DOUBLE, _rank_id_neighbours[d][(direction + 1) % 2], 99,
                      MPIDomain::sim_processor.comm, &recv_requests[d][direction]);
        }
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numrecv = numPartsToRecv[d][direction];
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);

            //将收到的电子云密度信息加到对应存储位置上
            pack::unpack_rho(d, direction, _atom->getAtomListRef(), recvbuf[direction], sendlist);
//            _atom->unpack_rho(d, direction, recvbuf[direction], sendlist);
            // 释放buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}

void Domain::sendDfEmbed(atom *_atom) {
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
            pack::pack_df(_atom->getAtomListRef(), sendbuf[direction],
                          _atom->inter_atom_list,
                          sendlist[iswap], intersendlist[iswap]);
//            _atom->pack_df(sendlist[iswap], intersendlist[iswap], sendbuf[direction]);
            iswap++;
        }

        // 与上下邻居通信
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numsend = numPartsToSend[d][direction];
            int numrecv;

            MPI_Isend(sendbuf[direction], numsend, MPI_DOUBLE, _rank_id_neighbours[d][direction], 99,
                      MPIDomain::sim_processor.comm, &send_requests[d][direction]);
            MPI_Probe(_rank_id_neighbours[d][(direction + 1) % 2], 99, MPIDomain::sim_processor.comm,
                      &status);//测试邻居是否有信息发送给本地
            MPI_Get_count(&status, MPI_DOUBLE, &numrecv);//得到要接收的粒子数目
            // 初始化接收缓冲区
            //依据得到发送方要发送粒子信息大小，初始化接收缓冲区
            recvbuf[direction] = new double[numrecv];
            numPartsToRecv[d][direction] = numrecv;
            MPI_Irecv(recvbuf[direction], numrecv, MPI_DOUBLE, _rank_id_neighbours[d][(direction + 1) % 2], 99,
                      MPIDomain::sim_processor.comm, &recv_requests[d][direction]);
        }

        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numrecv = numPartsToRecv[d][direction];
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);

            //将收到的嵌入能导数信息加到对应存储位置上
            pack::unpack_df(numrecv, _atom->getAtomListRef(), recvbuf[direction],
                            _atom->inter_atom_list,
                            recvlist[jswap], interrecvlist[jswap]);
//            _atom->unpack_df(numrecv, recvbuf[direction], recvlist[jswap], interrecvlist[jswap]);
            jswap++;

            // release memory of buffer
            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}

void Domain::sendForce(atom *_atom) {
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
            pack::pack_force(numPartsToSend[d][direction], _atom->getAtomListRef(),
                             sendbuf[direction], recvlist[iswap--]);
        }
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numsend = numPartsToSend[d][direction] * 3;
            int numrecv;

            MPI_Isend(sendbuf[direction], numsend, MPI_DOUBLE, _rank_id_neighbours[d][direction], 99,
                      MPIDomain::sim_processor.comm, &send_requests[d][direction]);
            MPI_Probe(_rank_id_neighbours[d][(direction + 1) % 2], 99, MPIDomain::sim_processor.comm,
                      &status);//测试邻居是否有信息发送给本地
            MPI_Get_count(&status, MPI_DOUBLE, &numrecv);//得到要接收的粒子数目
            // 初始化接收缓冲区
            //依据得到发送方要发送粒子信息大小，初始化接收缓冲区
            recvbuf[direction] = new double[numrecv];
            numPartsToRecv[d][direction] = numrecv;
            MPI_Irecv(recvbuf[direction], numrecv, MPI_DOUBLE, _rank_id_neighbours[d][(direction + 1) % 2], 99,
                      MPIDomain::sim_processor.comm, &recv_requests[d][direction]);
        }
        for (direction = LOWER; direction <= HIGHER; direction++) {
            int numrecv = numPartsToRecv[d][direction]; // todo remove not used variable.
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);

            //将收到的粒子位置信息加到对应存储位置上
            pack::unpack_force(d, direction, _atom->getAtomListRef(), recvbuf[direction], sendlist);

            delete[] sendbuf[direction];
            delete[] recvbuf[direction];
        }
    }
}
