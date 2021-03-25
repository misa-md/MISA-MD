//
// Created by genshen on 2021/3/24.
//

#ifndef MISA_MD_SEND_RECV_LISTS_H
#define MISA_MD_SEND_RECV_LISTS_H

#include <vector>
#include <comm/domain/domain.h>
#include <comm/domain/bcc_domain.h>

#include "../atom/atom_list.h"
#include "../atom/inter_atom_list.h"
#include "../lattice/lattice.h"
#include "../types/pre_define.h"

class SendRecvLists {
    friend class atom; // todo remove friend class
public:
    explicit SendRecvLists(AtomList &atom_list, InterAtomList &inter_atom_list) :
            atom_list(atom_list), inter_atom_list(inter_atom_list), intersendlist(6), interrecvlist(6) {}

    inline std::vector<std::vector<_type_atom_id> > &getSendList() { return sendlist; }

    inline std::vector<std::vector<_type_atom_id> > &getRecvList() { return recvlist; }

    inline _type_inter_buf &getInterSendList() { return intersendlist; }

    inline _type_inter_buf &getInterRecvList() { return interrecvlist; }

    void exchangeAtomFirst(comm::BccDomain *p_domain);

    void exchangeAtom(comm::BccDomain *p_domain);

    void exchangeInter(comm::Domain *p_domain);

    /**
     * setup ghost area for inter atoms.
     * send inter atoms in simulation area that are contributed to ghost area of other processes.
     * @param p_domain pointer of domain
     */
    void borderInter(comm::BccDomain *p_domain);

private:
    AtomList &atom_list;
    InterAtomList &inter_atom_list;

    // the array to record atoms that are out of box.
    std::vector<std::vector<_type_atom_id> > sendlist; // todo make it temp data
    std::vector<std::vector<_type_atom_id> > recvlist;

    _type_inter_buf intersendlist; // the atoms to be send to other processes as ghost.
    _type_inter_buf interrecvlist;
};


#endif //MISA_MD_SEND_RECV_LISTS_H
