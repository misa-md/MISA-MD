//
// Created by genshen on 2021/3/24.
//

#include <utils/mpi_domain.h>
#include <comm/comm.hpp>
#include <comm/preset/comm_forwarding_region.h>

#include "send_recv_lists.h"
#include "../pack/lat_particle_packer.h"
#include "../utils/mpi_data_types.h"


void SendRecvLists::exchangeAtomFirst(comm::BccDomain *p_domain) {
    sendlist.resize(6);
    recvlist.resize(6);
    for (unsigned short d = 0; d < DIMENSION; d++) {
        for (int direction = comm::DIR_LOWER; direction <= comm::DIR_HIGHER; direction++) {
            // 找到要发送给邻居的原子
            const int send_list_index = 2 * d + direction;
            std::vector<_type_atom_id> &dd_send_list = sendlist[send_list_index];
            comm::Region<comm::_type_lattice_size> region = comm::fwCommLocalRegion(p_domain, d, direction);
            for (int iz = region.z_low; iz < region.z_high; iz++) {
                for (int iy = region.y_low; iy < region.y_high; iy++) {
                    for (int ix = region.x_low; ix < region.x_high; ix++) {
                        dd_send_list.push_back(atom_list.lattice.IndexOf3DIndex(ix, iy, iz));
                    }
                }
            }
        }
    }

    LatPackerFirst lat_packer(*p_domain, atom_list, sendlist, recvlist);
    comm::neiSendReceive<LatParticleData>(&lat_packer,
                                          MPIDomain::toCommProcess(),
                                          mpi_types::_mpi_latParticle_data,
                                          p_domain->rank_id_neighbours);
}

void SendRecvLists::exchangeAtom(comm::BccDomain *p_domain) {
    LatPacker lat_packer(*p_domain, atom_list, sendlist, recvlist);
    comm::neiSendReceive<LatParticleData>(&lat_packer,
                                          MPIDomain::toCommProcess(),
                                          mpi_types::_mpi_latParticle_data,
                                          p_domain->rank_id_neighbours);
}
