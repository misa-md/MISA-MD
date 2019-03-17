//
// Created by genshen on 2019-03-16.
//

#include "df_embed_packer.h"
#include "pack.h"

DfEmbedPacker::DfEmbedPacker(AtomList &atom_list, InterAtomList &inter_atom_list,
                             std::vector<std::vector<_type_atom_id>> &send_list,
                             std::vector<std::vector<_type_atom_id>> &receive_list,
                             _type_inter_buf &inter_send_list,
                             _type_inter_buf &inter_receive_list)
        : atom_list(atom_list), inter_atom_list(inter_atom_list),
          send_list(send_list), receive_list(receive_list),
          inter_send_list(inter_send_list), inter_receive_list(inter_receive_list) {}

const unsigned long DfEmbedPacker::sendLength(const int dimension, const int direction) {
    const int index = 2 * dimension + direction;
    return send_list[index].size() + inter_send_list[index].size();;
}

void DfEmbedPacker::onSend(double *buffer, const unsigned long send_len,
                           const int dimension, const int direction) {
    const int index = 2 * dimension + direction;
    pack::pack_df(atom_list, buffer, &inter_atom_list,
                  send_list[index], inter_send_list[index]);
}

void DfEmbedPacker::onReceive(double *buffer, const unsigned long receive_len,
                              const int dimension, const int direction) {
    const int index = 2 * dimension + direction;
    //将收到的嵌入能导数信息加到对应存储位置上
    pack::unpack_df(receive_len, atom_list, buffer, &inter_atom_list,
                    receive_list[index], inter_receive_list[index]);
}
