//
// Created by genshen on 5/19/18.
//

#include <cmath>
#include <domain/domain.h>
#include <comm.hpp>
#include <pack/inter_border_packer.h>

#include "pack/inter_particle_packer.h"
#include "inter_atom_list.h"
#include "../utils/mpi_domain.h"
#include "../utils/mpi_data_types.h"

InterAtomList::InterAtomList() : nlocalinter(0), nghostinter(0),
                                 intersendlist(6), interrecvlist(6) {}

void InterAtomList::exchangeInter(comm::Domain *p_domain) {
    InterParticlePacker inter_packer(*p_domain, *this);
    comm::neiSendReceive<particledata>(&inter_packer,
                                       MPIDomain::toCommProcess(),
                                       mpi_types::_mpi_Particle_data,
                                       p_domain->rank_id_neighbours);
}

void InterAtomList::appendInter(_type_atom_id atom_id) {

}

void InterAtomList::borderInter(comm::Domain *p_domain) {
    InterBorderPacker border_packer(*p_domain, *this);
    comm::neiSendReceive<LatParticleData>(&border_packer,
                                          MPIDomain::toCommProcess(),
                                          mpi_types::_mpi_latParticle_data,
                                          p_domain->rank_id_neighbours);
}

void InterAtomList::addInterAtom(AtomElement &atom) {
    inter_list.push_back(atom);
    nlocalinter++;
}

void InterAtomList::addGhostAtom(AtomElement &ghost_atom) {
    inter_ghost_list.push_back(ghost_atom);
    nghostinter++;
}

_type_inter_list::iterator InterAtomList::removeInter(_type_inter_list::iterator inter_it) {
    // remove it.
    nlocalinter--;
    return inter_list.erase(inter_it);
}
