//
// Created by genshen on 2019-03-16.
//

#include <types_define.h>

#include "rho_packer.h"
#include "pack.h"

RhoPacker::RhoPacker(AtomList &atom_list, std::vector<std::vector<_type_atom_id>> &send_list,
                     std::vector<std::vector<_type_atom_id>> &receive_list)
        : atom_list(atom_list), send_list(send_list), receive_list(receive_list) {}

const unsigned long RhoPacker::sendLength(const int dimension, const int direction) {
    // dimension: 2   1    0
    // direction:L H  L H  L H
    // index:    5 4  3 2  1 0
    const int index = 2 * dimension + (direction == LOWER ? 1 : 0);
    return receive_list[index].size();
}

void RhoPacker::onSend(double *buffer, const unsigned long send_len,
                       const int dimension, const int direction) {
    const int index = 2 * dimension + (direction == LOWER ? 1 : 0);
    pack::pack_rho(send_len, atom_list,
                   buffer, receive_list[index]);
}

void RhoPacker::onReceive(double buffer[], const unsigned long receive_len,
                          const int dimension, const int direction) {
    //将收到的电子云密度信息加到对应存储位置上
    pack::unpack_rho(dimension, direction, atom_list, buffer, send_list);
}
