//
// Created by genshen on 2019-03-17.
//

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
    getIntertosend(&domain, dimension, direction,
                   domain.meas_ghost_length[dimension],
                   inter_atom_list.intersendlist[index]);
    return inter_atom_list.intersendlist[index].size();
}

void InterBorderPacker::onSend(LatParticleData buffer[],
                               const unsigned long send_len,
                               const int dimension,
                               const int direction) {
    double shift = 0.0;
    if (domain.grid_coord_sub_box[dimension] == 0 && direction == LOWER) {
        shift = domain.meas_global_length[dimension];
    }
    if (domain.grid_coord_sub_box[dimension] == domain.grid_size[dimension] - 1 && direction == HIGHER) {
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
    for (int i = 0; i < receive_len; i++) {
        ele.type = buffer[i].type; // todo check type invalid
        ele.x[0] = buffer[i].r[0];
        ele.x[1] = buffer[i].r[1];
        ele.x[2] = buffer[i].r[2];
        if (ele.x[0] >= domain.meas_ghost_region.low[0] &&
            ele.x[0] < domain.meas_ghost_region.high[0] &&
            ele.x[1] >= domain.meas_ghost_region.low[1] &&
            ele.x[1] < domain.meas_ghost_region.high[1] &&
            ele.x[2] >= domain.meas_ghost_region.low[2] &&
            ele.x[2] < domain.meas_ghost_region.high[2]) {
            inter_atom_list.inter_ghost_list.push_back(ele);
            inter_atom_list.nghostinter++;
        } else {
            // todo warning
            inter_atom_list.interrecvlist[index][i] = nullptr;
        }
    }
}

void InterBorderPacker::getIntertosend(const comm::Domain *p_domain, int d, int direction,
                                       double ghostlengh, std::vector<AtomElement *> &sendlist) {
    double low, high;
    if (d == 0) {
        if (direction == 0) {
            low = p_domain->meas_sub_box_region.x_low;
            high = p_domain->meas_sub_box_region.x_low + ghostlengh;
            for (AtomElement &inter_ref : inter_atom_list.inter_list) {
                if (inter_ref.x[0] < high && inter_ref.x[0] >= low) {
                    sendlist.push_back(&inter_ref);
                }
            }
            for (AtomElement &ghost_ref :inter_atom_list.inter_ghost_list) {
                if (ghost_ref.x[0] < high && ghost_ref.x[0] >= low) {
                    sendlist.push_back(&ghost_ref);
                }
            }
        } else {
            low = p_domain->meas_sub_box_region.x_high - ghostlengh;
            high = p_domain->meas_sub_box_region.x_high;
            for (AtomElement &inter_ref :inter_atom_list.inter_list) {
                if (inter_ref.x[0] <= high && inter_ref.x[0] > low) {
                    sendlist.push_back(&inter_ref);
                }
            }
            for (AtomElement &ghost_ref :inter_atom_list.inter_ghost_list) {
                if (ghost_ref.x[0] <= high && ghost_ref.x[0] > low) {
                    sendlist.push_back(&ghost_ref);
                }
            }
        }
    } else if (d == 1) {
        if (direction == 0) {
            low = p_domain->meas_sub_box_region.y_low;
            high = p_domain->meas_sub_box_region.y_low + ghostlengh;
            for (AtomElement &inter_ref :inter_atom_list.inter_list) {
                if (inter_ref.x[1] < high && inter_ref.x[1] >= low) {
                    sendlist.push_back(&inter_ref);
                }
            }
            for (AtomElement &ghost_ref :inter_atom_list.inter_ghost_list) {
                if (ghost_ref.x[1] < high && ghost_ref.x[1] >= low) {
                    sendlist.push_back(&ghost_ref);
                }
            }
        } else {
            low = p_domain->meas_sub_box_region.y_high - ghostlengh;
            high = p_domain->meas_sub_box_region.y_high;
            for (AtomElement &inter_ref :inter_atom_list.inter_list) {
                if (inter_ref.x[1] <= high && inter_ref.x[1] > low) {
                    sendlist.push_back(&inter_ref);
                }
            }
            for (AtomElement &ghost_ref :inter_atom_list.inter_ghost_list) {
                if (ghost_ref.x[1] <= high && ghost_ref.x[1] > low) {
                    sendlist.push_back(&ghost_ref);
                }
            }
        }
    } else {
        if (direction == 0) {
            low = p_domain->meas_sub_box_region.z_low;
            high = p_domain->meas_sub_box_region.z_low + ghostlengh;
            for (AtomElement &inter_ref :inter_atom_list.inter_list) {
                if (inter_ref.x[2] < high && inter_ref.x[2] >= low) {
                    sendlist.push_back(&inter_ref);
                }
            }
            for (AtomElement &ghost_ref :inter_atom_list.inter_ghost_list) {
                if (ghost_ref.x[2] < high && ghost_ref.x[2] >= low) {
                    sendlist.push_back(&ghost_ref);
                }
            }
        } else {
            low = p_domain->meas_sub_box_region.z_high - ghostlengh;
            high = p_domain->meas_sub_box_region.z_high;
            for (AtomElement &inter_ref :inter_atom_list.inter_list) {
                if (inter_ref.x[2] <= high && inter_ref.x[2] > low) {
                    sendlist.push_back(&inter_ref);
                }
            }
            for (AtomElement &ghost_ref :inter_atom_list.inter_ghost_list) {
                if (ghost_ref.x[2] <= high && ghost_ref.x[2] > low) {
                    sendlist.push_back(&ghost_ref);
                }
            }
        }
    }
}