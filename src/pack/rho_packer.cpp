//
// Created by genshen on 2019-03-16.
//

#include <comm/types_define.h>

#include "rho_packer.h"

RhoPacker::RhoPacker(AtomList &atom_list, std::vector<std::vector<_type_atom_id>> &send_list,
                     std::vector<std::vector<_type_atom_id>> &receive_list)
        : atom_list(atom_list), send_list(send_list), receive_list(receive_list) {}

const unsigned long RhoPacker::sendLength(const int dimension, const int direction) {
    // dimension: 2   1    0
    // direction:L H  L H  L H
    // index:    5 4  3 2  1 0
    const int index = 2 * dimension + (direction == comm::DIR_LOWER ? 1 : 0);
    return receive_list[index].size();
}

void RhoPacker::onSend(double *buffer, const unsigned long send_len,
                       const int dimension, const int direction) {
    const int index = 2 * dimension + (direction == comm::DIR_LOWER ? 1 : 0);
    std::vector<_type_atom_id> &recvlist = receive_list[index];
    int j, m = 0;
    for (int i = 0; i < send_len; i++) {
        j = recvlist[i];
        AtomElement &atom = atom_list._atoms.getAtomEleByLinearIndex(j);
        buffer[m++] = atom.rho;
    }
}

void RhoPacker::onReceive(double buffer[], const unsigned long receive_len,
                          const int dimension, const int direction) {
    //将收到的电子云密度信息加到对应存储位置上
    double *buf = buffer;
    // flip the direction
    const int list_index = 2 * dimension + (direction == comm::DIR_LOWER ? comm::DIR_HIGHER : comm::DIR_LOWER);

    int j, m = 0;
    for (int i = 0; i < send_list[list_index].size(); i++) {
        j = send_list[list_index][i];
        AtomElement &atom_ = atom_list._atoms.getAtomEleByLinearIndex(j);
        atom_.rho += buf[m++];
    }
}
