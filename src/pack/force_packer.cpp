//
// Created by genshen on 2019-03-16.
//

#include "force_packer.h"

ForcePacker::ForcePacker(AtomList &atom_list, std::vector<std::vector<_type_atom_id>> &send_list,
                         std::vector<std::vector<_type_atom_id>> &receive_list)
        : atom_list(atom_list), send_list(send_list), receive_list(receive_list) {}

const unsigned long ForcePacker::sendLength(const int dimension, const int direction) {
    const int index = 2 * dimension + (direction == comm::DIR_LOWER ? 1 : 0);
    return receive_list[index].size() * 3;
}

void ForcePacker::onSend(double *buffer, const unsigned long send_len, const int dimension, const int direction) {
    const int index = 2 * dimension + (direction == comm::DIR_LOWER ? 1 : 0); // fixme possible bug
    std::vector<_type_atom_id> &recvlist = receive_list[index];
    int j, m = 0;
    for (int i = 0; i < send_len / 3; i++) {
        j = recvlist[i];
        AtomElement &atom = atom_list.getAtomEleByLinearIndex(j);
        buffer[m++] = atom.f[0];
        buffer[m++] = atom.f[1];
        buffer[m++] = atom.f[2];
    }
}

void ForcePacker::onReceive(double *buffer, const unsigned long receive_len, const int dimension, const int direction) {
    //将收到的粒子位置信息加到对应存储位置上
    const int list_index = 2 * dimension + (direction == comm::DIR_LOWER ? comm::DIR_HIGHER : comm::DIR_LOWER); // Flip the direction
    int j, m = 0;
    for (int i = 0; i < send_list[list_index].size(); i++) {
        j = send_list[list_index][i];
        AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
        atom_.f[0] += buffer[m++];
        atom_.f[1] += buffer[m++];
        atom_.f[2] += buffer[m++];
    }
}
