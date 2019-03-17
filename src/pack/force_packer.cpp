//
// Created by genshen on 2019-03-16.
//

#include "force_packer.h"
#include "pack.h"

ForcePacker::ForcePacker(AtomList &atom_list, std::vector<std::vector<_type_atom_id>> &send_list,
                         std::vector<std::vector<_type_atom_id>> &receive_list)
        : atom_list(atom_list), send_list(send_list), receive_list(receive_list) {}

const unsigned long ForcePacker::sendLength(const int dimension, const int direction) {
    const int index = 2 * dimension + (direction == LOWER ? 1 : 0);
    return receive_list[index].size() * 3;
}

void ForcePacker::onSend(double *buffer, const unsigned long send_len, const int dimension, const int direction) {
    const int index = 2 * dimension + (direction == LOWER ? 1 : 0);
    pack::pack_force(send_len / 3, atom_list, buffer, receive_list[index]);
}

void ForcePacker::onReceive(double *buffer, const unsigned long receive_len, const int dimension, const int direction) {
//    const int index = 2 * dimension + (direction == LOWER ? 1 : 0);
    //将收到的粒子位置信息加到对应存储位置上
    pack::unpack_force(dimension, direction, atom_list, buffer, send_list);
}
