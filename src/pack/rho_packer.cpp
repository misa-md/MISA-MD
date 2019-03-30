//
// Created by genshen on 2019-03-16.
//

#include <types_define.h>

#include "rho_packer.h"

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
    std::vector<_type_atom_id> &recvlist = receive_list[index];
    int j, m = 0;
    for (int i = 0; i < send_len; i++) {
        j = recvlist[i];
        AtomElement &atom = atom_list.getAtomEleByLinearIndex(j);
        buffer[m++] = atom.rho;
    }
}

void RhoPacker::onReceive(double buffer[], const unsigned long receive_len,
                          const int dimension, const int direction) {
    //将收到的电子云密度信息加到对应存储位置上
    const int d = dimension;
    double *buf = buffer;
    std::vector<std::vector<_type_atom_id> > &sendlist = send_list;
    int j, m = 0;
    if (d == 0) {
        if (direction == 0) {
            for (int i = 0; i < sendlist[1].size(); i++) {
                j = sendlist[1][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[0].size(); i++) {
                j = sendlist[0][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        }
    } else if (d == 1) {
        if (direction == 0) {
            for (int i = 0; i < sendlist[3].size(); i++) {
                j = sendlist[3][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[2].size(); i++) {
                j = sendlist[2][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        }
    } else {
        if (direction == 0) {
            for (int i = 0; i < sendlist[5].size(); i++) {
                j = sendlist[5][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[4].size(); i++) {
                j = sendlist[4][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        }
    }
}
