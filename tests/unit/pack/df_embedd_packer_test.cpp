//
// Created by genshen on 2021/2/18.
//

#include <gtest/gtest.h>
#include <comm/comm.hpp>
#include <utils/mpi_domain.h>
#include <comm/preset/comm_forwarding_region.h>
#include "pack/df_embed_packer.h"

#include "../fixtures/world_builder_test_fixture.h"

//@MPI; @ranks=4,6
TEST_F(WorldBuilderTestFixture, df_embedd_packer_test) {
    // set inter sending and receiving list to nothing.
    _atom->p_send_recv_list->getInterSendList().clear();
    _atom->p_send_recv_list->getInterRecvList().clear();
    _atom->p_send_recv_list->getInterSendList().resize(6);
    _atom->p_send_recv_list->getInterRecvList().resize(6);

    std::vector<std::vector<_type_atom_id> > &send_list = _atom->p_send_recv_list->getSendList();
    std::vector<std::vector<_type_atom_id> > &recv_list = _atom->p_send_recv_list->getRecvList();
    send_list.clear();
    recv_list.clear();
    send_list.resize(6);
    recv_list.resize(6);

    MPIDomain::sim_processor = kiwi::mpiUtils::global_process;

    // build send list and receive list
    for (unsigned short d = 0; d < DIMENSION; d++) {
        for (int direction = comm::DIR_LOWER; direction <= comm::DIR_HIGHER; direction++) {
            std::vector<_type_atom_id> &dd_send_list = send_list[2 * d + direction];
            comm::Region<comm::_type_lattice_size> region = comm::fwCommLocalSendRegion(
                    p_domain->dbx_lattice_size_ghost,
                    p_domain->dbx_local_sub_box_lattice_region,
                    d, direction);
            for (int iz = region.z_low; iz < region.z_high; iz++) {
                for (int iy = region.y_low; iy < region.y_high; iy++) {
                    for (int ix = region.x_low; ix < region.x_high; ix++) {
                        dd_send_list.push_back(_atom->getAtomList()->lattice.IndexOf3DIndex(ix, iy, iz));
                    }
                }
            }
        }
    }
    for (unsigned short d = 0; d < DIMENSION; d++) {
        for (int direction = comm::DIR_LOWER; direction <= comm::DIR_HIGHER; direction++) {
            std::vector<_type_atom_id> &dd_recv_list = recv_list[2 * d + direction];
            comm::Region<comm::_type_lattice_size> region = comm::fwCommLocalRecvRegion(
                    p_domain->dbx_lattice_size_ghost,
                    p_domain->dbx_local_sub_box_lattice_region,
                    d, direction);
            for (int iz = region.z_low; iz < region.z_high; iz++) {
                for (int iy = region.y_low; iy < region.y_high; iy++) {
                    for (int ix = region.x_low; ix < region.x_high; ix++) {
                        dd_recv_list.push_back(_atom->getAtomList()->lattice.IndexOf3DIndex(ix, iy, iz));
                    }
                }
            }
        }
    }

    // set df values
    _atom->atom_list->foreachSubBoxAtom([&](const _type_atom_index gid) {
        MD_LOAD_ATOM_VAR(element, _atom->atom_list, gid);
        MD_SET_ATOM_DF(element, gid, MPIDomain::sim_processor.own_rank);
    });

    DfEmbedPacker packer(_atom->getAtomListRef(),
                         _atom->p_send_recv_list->getSendList(), _atom->p_send_recv_list->getRecvList(),
                         _atom->p_send_recv_list->getInterSendList(),
                         _atom->p_send_recv_list->getInterRecvList());

    comm::neiSendReceive<double>(&packer, MPIDomain::toCommProcess(), MPI_DOUBLE, p_domain->rank_id_neighbours);

    const _type_atom_index gid1 = _atom->atom_list->_atoms.getAtomIndex(
            0,
            p_domain->dbx_local_sub_box_lattice_region.y_high / 2,
            p_domain->dbx_local_sub_box_lattice_region.z_high / 2);
    MD_LOAD_ATOM_VAR(_atom_1, _atom->atom_list, gid1);
    double df1 = MD_GET_ATOM_DF(_atom_1, gid1);
    EXPECT_EQ(df1, p_domain->rank_id_neighbours[0][0]);

    const _type_atom_index gid2 = _atom->atom_list->_atoms.getAtomIndex(
            p_domain->dbx_local_sub_box_lattice_region.x_high + 1,
            p_domain->dbx_local_sub_box_lattice_region.y_high / 2,
            p_domain->dbx_local_sub_box_lattice_region.z_high / 2);
    MD_LOAD_ATOM_VAR(_atom_2, _atom->atom_list, gid2);
    double df2 = MD_GET_ATOM_DF(_atom_2, gid2);
    EXPECT_EQ(df2, p_domain->rank_id_neighbours[0][1]);
}
