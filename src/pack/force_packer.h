//
// Created by genshen on 2019-03-16.
//

#ifndef CRYSTAL_MD_FORCE_PACKER_H
#define CRYSTAL_MD_FORCE_PACKER_H

#include <vector>
#include <mpi.h>
#include <comm/packer.h>
#include "atom/atom_list.h"

/**
 * @brief force data packer for communicating with neighbour processes.
 */
class ForcePacker : public comm::Packer<double> {
public:
    ForcePacker(AtomList &atom_list,
                std::vector<std::vector<_type_atom_id>> &send_list,
                std::vector<std::vector<_type_atom_id>> &receive_list);

    const unsigned long sendLength(const int dimension, const int direction) override;

    void onSend(double buffer[], const unsigned long send_len,
                const int dimension, const int direction) override;

    void onReceive(double buffer[], const unsigned long receive_len,
                   const int dimension, const int direction) override;

private:
    AtomList &atom_list;
    std::vector<std::vector<_type_atom_id>> &send_list;
    std::vector<std::vector<_type_atom_id>> &receive_list;
};


#endif //CRYSTAL_MD_FORCE_PACKER_H
