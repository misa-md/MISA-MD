//
// Created by genshen on 2019-03-16.
//

#ifndef MISA_MD_DF_EMBED_PACKER_H
#define MISA_MD_DF_EMBED_PACKER_H


#include <vector>
#include <comm/packer.h>
#include "atom/inter_atom_list.h"
#include "atom/atom_list.h"

/**
 * @brief df embed data packer for communicating with neighbour processes.
 */
class DfEmbedPacker : public comm::Packer<double> {
public:
    DfEmbedPacker(AtomList &atom_list,
                  std::vector<std::vector<_type_atom_id>> &send_list,
                  std::vector<std::vector<_type_atom_id>> &receive_list,
                  _type_inter_buf &inter_send_list,
                  _type_inter_buf &inter_receive_list);

    const unsigned long sendLength(const int dimension, const int direction) override;

    void onSend(double buffer[], const unsigned long send_len,
                const int dimension, const int direction) override;

    void onReceive(double buffer[], const unsigned long receive_len,
                   const int dimension, const int direction) override;

private:
    AtomList &atom_list;
    std::vector<std::vector<_type_atom_id>> &send_list;
    std::vector<std::vector<_type_atom_id>> &receive_list;
    _type_inter_buf &inter_send_list;
    _type_inter_buf &inter_receive_list;
};


#endif //MISA_MD_DF_EMBED_PACKER_H
