//
// Created by genshen on 2019-03-17.
//

#include <cstring>
#include <logs/logs.h>
#include <comm_forwarding_region.h>
#include <atom/ws_utils.h>
#include "inter_border_packer.h"


InterBorderPacker::InterBorderPacker(const comm::Domain &domain, InterAtomList &inter_atom_list)
        : domain(domain), inter_atom_list(inter_atom_list) {
    inter_atom_list.intersendlist.clear();
    inter_atom_list.interrecvlist.clear();
    inter_atom_list.intersendlist.resize(6);
    inter_atom_list.interrecvlist.resize(6);
}

const unsigned long InterBorderPacker::sendLength(const int dimension, const int direction) {
    const int index = 2 * dimension + direction;
    std::vector<AtomElement *> &sendlist = inter_atom_list.intersendlist[index];
    // before x dimension communication, ghost list is empty.
    comm::Region<comm::_type_lattice_size> region = comm::fwCommLocalRegion(&domain, dimension, direction);
    _type_atom_index coords[DIMENSION] = {0, 0, 0};

    for (AtomElement &inter_ref : inter_atom_list.inter_list) {
        // get the lattice coordinate the inter atom belongs to.
        ws::getNearLatCoord(inter_ref, &domain, coords);
        if (region.isIn(coords[0], coords[1], coords[2])) {
            sendlist.push_back(&inter_ref);
        }
    }
    for (AtomElement &ghost_ref :inter_atom_list.inter_ghost_list) {
        ws::getNearLatCoord(ghost_ref, &domain, coords);
        if (region.isIn(coords[0], coords[1], coords[2])) {
            sendlist.push_back(&ghost_ref);
        }
    }
    return sendlist.size();
}

void InterBorderPacker::onSend(LatParticleData buffer[],
                               const unsigned long send_len,
                               const int dimension,
                               const int direction) {
    double shift = 0.0;
    if (domain.grid_coord_sub_box[dimension] == 0 && direction == comm::DIR_LOWER) {
        shift = domain.meas_global_length[dimension];
    }
    if (domain.grid_coord_sub_box[dimension] == domain.grid_size[dimension] - 1 && direction == comm::DIR_HIGHER) {
        shift = -((domain.meas_global_length[dimension]));
    }

    // fixme pack atom id?
    const int index = 2 * dimension + direction;
    switch (dimension) {
        case 0:
            for (int i = 0; i < send_len; i++) {
                buffer[i].type = inter_atom_list.intersendlist[index][i]->type;
                buffer[i].r[0] = inter_atom_list.intersendlist[index][i]->x[0] + shift;
                buffer[i].r[1] = inter_atom_list.intersendlist[index][i]->x[1];
                buffer[i].r[2] = inter_atom_list.intersendlist[index][i]->x[2];
            }
            break;
        case 1:
            for (int i = 0; i < send_len; i++) {
                buffer[i].type = inter_atom_list.intersendlist[index][i]->type;
                buffer[i].r[0] = inter_atom_list.intersendlist[index][i]->x[0];
                buffer[i].r[1] = inter_atom_list.intersendlist[index][i]->x[1] + shift;
                buffer[i].r[2] = inter_atom_list.intersendlist[index][i]->x[2];
            }
            break;
        case 2:
            for (int i = 0; i < send_len; i++) {
                buffer[i].type = inter_atom_list.intersendlist[index][i]->type;
                buffer[i].r[0] = inter_atom_list.intersendlist[index][i]->x[0];
                buffer[i].r[1] = inter_atom_list.intersendlist[index][i]->x[1];
                buffer[i].r[2] = inter_atom_list.intersendlist[index][i]->x[2] + shift;
            }
            break;
        default:
            break;
    }
}

void InterBorderPacker::onReceive(LatParticleData buffer[],
                                  const unsigned long receive_len,
                                  const int dimension,
                                  const int direction) {
    const int index = 2 * dimension + direction;
    inter_atom_list.interrecvlist[index].resize(receive_len);

    AtomElement ele;
    memset(&ele, 0, sizeof(AtomElement)); // set f,v,rho,df to 0.
    for (int i = 0; i < receive_len; i++) {
        // id is not necessary for ghost inter atoms.
        ele.type = buffer[i].type; // all type are valid(inter atoms)
        ele.x[0] = buffer[i].r[0];
        ele.x[1] = buffer[i].r[1];
        ele.x[2] = buffer[i].r[2];
        // todo: check if the atom is in ghost area using lattice coord.
        // kiwi::logs::w("inter", "unexpected atom.\n");
        inter_atom_list.interrecvlist[index][i] = inter_atom_list.addGhostAtom(ele);
    }
}
