//
// Created by genshen on 2019-03-16.
//

#ifndef CRYSTAL_MD_RHO_PACKER_H
#define CRYSTAL_MD_RHO_PACKER_H

#include <comm/packer.h>
#include "atom/atom_list.h"

/**
 * @brief rho data packer for communicating with neighbour processes.
 */
class RhoPacker : public comm::Packer<double> {
public:
    RhoPacker(AtomList &atom_list,
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


#endif //CRYSTAL_MD_RHO_PACKER_H
