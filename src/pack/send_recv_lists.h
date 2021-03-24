//
// Created by genshen on 2021/3/24.
//

#ifndef MISA_MD_SEND_RECV_LISTS_H
#define MISA_MD_SEND_RECV_LISTS_H

#include <vector>
#include <comm/domain/domain.h>
#include <comm/domain/bcc_domain.h>

#include "../atom/atom_list.h"
#include "../lattice/lattice.h"
#include "../types/pre_define.h"

class SendRecvLists {
    friend class atom; // todo remove friend class
public:
    explicit SendRecvLists(AtomList &atom_list) : atom_list(atom_list) {}

    inline std::vector<std::vector<_type_atom_id> > &getSendList() { return sendlist; }

    inline std::vector<std::vector<_type_atom_id> > &getRecvList() { return recvlist; }

    void exchangeAtomFirst(comm::BccDomain *p_domain);

    void exchangeAtom(comm::BccDomain *p_domain);

private:
    AtomList &atom_list;

    // the array to record atoms that are out of box.
    std::vector<std::vector<_type_atom_id> > sendlist; // todo make it temp data
    std::vector<std::vector<_type_atom_id> > recvlist;
};


#endif //MISA_MD_SEND_RECV_LISTS_H
