//
// Created by genshen on 2018-12-28.
// Updated by genshen on 2019-03-17 to add libcomm.
//

#include <cstring>
#include <cassert>
#include <lattice/ws_utils.h>
#include <logs/logs.h>
#include "atom/inter_atom_list.h"
#include "inter_particle_packer.h"


InterParticlePacker::InterParticlePacker(const comm::Domain &domain, InterAtomList &inter_atom_list)
        : domain(domain), inter_atom_list(inter_atom_list) {}

const unsigned long InterParticlePacker::sendLength(const int dimension, const int direction) {
    const static box::_type_flag_32 out_box_flags[DIMENSION][2] = {
            {box::OUT_BOX_X_LITTER, box::OUT_BOX_X_BIG},
            {box::OUT_BOX_Y_LITTER, box::OUT_BOX_Y_BIG},
            {box::OUT_BOX_Z_LITTER, box::OUT_BOX_Z_BIG},
    };
    unsigned int num_to_send = 0;
    for (_type_inter_list::iterator inter_it = inter_atom_list.inter_list.begin();
         inter_it != inter_atom_list.inter_list.end(); inter_it++) {
#ifdef MD_DEV_MODE
        {
            const box::_type_flag_32 flag = ws::isOutBox(*inter_it, &domain);
            if (dimension == 1) { // y dimension
                // In y dimension communication, atoms in inter atoms list
                // can not be out of x boundary of simulation box.
                // Because we have a communication in x dimension before.
                assert((flag & box::OUT_BOX_X_LITTER) == 0);
                assert((flag & box::OUT_BOX_X_BIG) == 0);
            }
            if (dimension == 2) { // z dimension
                // In y dimension communication, atoms in inter atoms list
                // can not be out of x and y boundary of simulation box due to preview communications.
                assert((flag & box::OUT_BOX_X_LITTER) == 0);
                assert((flag & box::OUT_BOX_X_BIG) == 0);
                assert((flag & box::OUT_BOX_Y_LITTER) == 0);
                assert((flag & box::OUT_BOX_Y_BIG) == 0);
            }
        }
#endif
        // we assume that, a atom cannot cross 2 or more than 2 sub-boxes
        if (ws::isOutBox(*inter_it, &domain) & out_box_flags[dimension][direction]) {
            num_to_send++;
        }
    }
    return num_to_send;
}

void InterParticlePacker::onSend(particledata *buffer, const unsigned long send_len,
                                 const int dimension, const int direction) {
    const static box::_type_flag_32 out_box_flags[DIMENSION][2] = {
            {box::OUT_BOX_X_LITTER, box::OUT_BOX_X_BIG},
            {box::OUT_BOX_Y_LITTER, box::OUT_BOX_Y_BIG},
            {box::OUT_BOX_Z_LITTER, box::OUT_BOX_Z_BIG},
    };

    double offset[DIMENSION] = {0.0};
    // periodic boundary
    if (domain.grid_coord_sub_box[dimension] == 0 && direction == comm::DIR_LOWER) {
        offset[dimension] = domain.meas_global_length[dimension];
    }
    if (domain.grid_coord_sub_box[dimension] == domain.grid_size[dimension] - 1 && direction == comm::DIR_HIGHER) {
        offset[dimension] = -((domain.meas_global_length[dimension]));
    }

    _type_inter_list &inter_list = inter_atom_list.inter_list;
    unsigned long i = 0; // todo type
    for (_type_inter_list::iterator inter_it = inter_list.begin(); inter_it != inter_list.end();) {
        // we assume that, a atom cannot cross 2 or more than 2 sub-boxes
        if (ws::isOutBox(*inter_it, &domain) & out_box_flags[dimension][direction]) {
            buffer[i].id = inter_it->id;
            buffer[i].type = inter_it->type;
            buffer[i].r[0] = inter_it->x[0] + offset[0];
            buffer[i].r[1] = inter_it->x[1] + offset[1];
            buffer[i].r[2] = inter_it->x[2] + offset[2];
            buffer[i].v[0] = inter_it->v[0];
            buffer[i].v[1] = inter_it->v[1];
            buffer[i].v[2] = inter_it->v[2];
            // remove the inter atom.
            // exchange the atom at end of vector to the position of atom (inter[j]) te be removed .
            inter_it = inter_atom_list.removeInter(inter_it);
            i++;
        } else {
            inter_it++;
        }
    }
}

void InterParticlePacker::onReceive(particledata *buffer, const unsigned long receive_len,
                                    const int dimension, const int direction) {
    AtomElement atom{};
    memset(&atom, 0, sizeof(AtomElement)); // set f,rho,df to 0.
    for (int i = 0; i < receive_len; i++) {
        atom.id = buffer[i].id;
        atom.type = buffer[i].type;
        atom.x[0] = buffer[i].r[0];
        atom.x[1] = buffer[i].r[1];
        atom.x[2] = buffer[i].r[2];
        atom.v[0] = buffer[i].v[0];
        atom.v[1] = buffer[i].v[1];
        atom.v[2] = buffer[i].v[2];
        inter_atom_list.addInterAtom(atom);
    }
}
