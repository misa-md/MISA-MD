//
// Created by genshen on 2019-03-16.
//

#include <mpi.h>
#include "df_embed_packer.h"

DfEmbedPacker::DfEmbedPacker(AtomList &atom_list,
                             std::vector<std::vector<_type_atom_id>> &send_list,
                             std::vector<std::vector<_type_atom_id>> &receive_list,
                             _type_inter_buf &inter_send_list,
                             _type_inter_buf &inter_receive_list)
        : atom_list(atom_list), send_list(send_list), receive_list(receive_list),
          inter_send_list(inter_send_list), inter_receive_list(inter_receive_list) {}

const unsigned long DfEmbedPacker::sendLength(const int dimension, const int direction) {
    const int index = 2 * dimension + direction;
    return send_list[index].size() + inter_send_list[index].size();;
}

void DfEmbedPacker::onSend(double *buffer, const unsigned long send_len,
                           const int dimension, const int direction) {
    const int index = 2 * dimension + direction;
    std::vector<_type_atom_id> &sendlist = send_list[index];
    std::vector<AtomElement *> &intersendlist = inter_send_list[index];
    int j, m = 0;
    int n = sendlist.size();
    for (int i = 0; i < n; i++) {
        j = sendlist[i];
        AtomElement &atom = atom_list._atoms.getAtomEleByLinearIndex(j);
        buffer[m++] = atom.df;
    }
    n = intersendlist.size();
    for (int i = 0; i < n; i++) {
        //   j = intersendlist[i];
        buffer[m++] = intersendlist[i]->df;
    }
}

void DfEmbedPacker::onReceive(double *buffer, const unsigned long receive_len,
                              const int dimension, const int direction) {
    const int index = 2 * dimension + direction;
    //将收到的嵌入能导数信息加到对应存储位置上
    std::vector<_type_atom_id> &recvlist = receive_list[index];
    std::vector<AtomElement *> &interrecvlist = inter_receive_list[index];

    long kk;
    auto m = 0;
    if (receive_len != (recvlist.size() + interrecvlist.size())) {
        printf("wrong number of dfembed recv!!!");
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
    unsigned long len1 = recvlist.size();
    for (int i = 0; i < len1; i++) {
        kk = recvlist[i];
        AtomElement &atom_ = atom_list._atoms.getAtomEleByLinearIndex(kk);
        atom_.df = buffer[m++];
    }
    unsigned long len2 = interrecvlist.size();
    for (int i = 0; i < len2; i++) {
        if (interrecvlist[i] == nullptr) {
            m++; // todo error
            continue;
        }
        interrecvlist[i]->df = buffer[m++];
    }
}
