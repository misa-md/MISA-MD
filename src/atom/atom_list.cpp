//
// Created by genshen on 5/8/18.
//

#include <utils/mpi_domain.h>
#include <pack/lat_particle_packer.h>
#include <comm.hpp>
#include <comm_forwarding_region.h>
#include "atom_list.h"
#include "../utils/mpi_data_types.h"

AtomList::AtomList(_type_atom_count size_x, _type_atom_count size_y, _type_atom_count size_z,
                   _type_atom_count size_sub_box_x, _type_atom_count size_sub_box_y, _type_atom_count size_sub_box_z,
                   _type_atom_count ghost_count_x, _type_atom_count ghost_count_y, _type_atom_count ghost_count_z) :
        lattice(size_x, size_y, size_z, size_sub_box_x, size_sub_box_y, size_sub_box_z,
                ghost_count_x, ghost_count_y, ghost_count_z) {

    _atoms = new AtomElement **[size_z];
    for (_type_atom_count z = 0; z < size_z; z++) {
        _atoms[z] = new AtomElement *[size_y];
        for (_type_atom_count y = 0; y < size_y; y++) {
            _atoms[z][y] = new AtomElement[size_x];
        }
    }
}

AtomList::~AtomList() {
    for (_type_atom_count z = 0; z < lattice._size_z; z++) {
        for (_type_atom_count y = 0; y < lattice._size_y; y++) {
            delete[] _atoms[z][y];
        }
        delete[] _atoms[z];
    }
    delete[] _atoms;
}

void AtomList::exchangeAtomFirst(comm::BccDomain *p_domain) {
    sendlist.resize(6);
    recvlist.resize(6);
    for (unsigned short d = 0; d < DIMENSION; d++) {
        for (int direction = comm::DIR_LOWER; direction <= comm::DIR_HIGHER; direction++) {
            // 找到要发送给邻居的原子
            const int send_list_index = 2 * d + direction;
            std::vector<_type_atom_id> &dd_send_list = sendlist[send_list_index];
            comm::Region<comm::_type_lattice_size> region = comm::fwCommLocalRegion(p_domain, d, direction);
            for (int iz = region.z_low; iz < region.z_high; iz++) {
                for (int iy = region.y_low; iy < region.y_high; iy++) {
                    for (int ix = region.x_low; ix < region.x_high; ix++) {
                        dd_send_list.push_back(lattice.IndexOf3DIndex(ix, iy, iz));
                    }
                }
            }
        }
    }

    LatPackerFirst lat_packer(*p_domain, *this, sendlist, recvlist);
    comm::neiSendReceive<LatParticleData>(&lat_packer,
                                          MPIDomain::toCommProcess(),
                                          mpi_types::_mpi_latParticle_data,
                                          p_domain->rank_id_neighbours);
}

void AtomList::exchangeAtom(comm::BccDomain *p_domain) {
    LatPacker lat_packer(*p_domain, *this, sendlist, recvlist);
    comm::neiSendReceive<LatParticleData>(&lat_packer,
                                          MPIDomain::toCommProcess(),
                                          mpi_types::_mpi_latParticle_data,
                                          p_domain->rank_id_neighbours);
}

bool AtomList::isBadList(comm::Domain domain) {
    const double delta = 10;
    for (long z = lattice.purge_ghost_count_z; z < lattice._size_sub_box_z + lattice.purge_ghost_count_z; z++) {
        for (long y = lattice.purge_ghost_count_y; y < lattice._size_sub_box_y + lattice.purge_ghost_count_y; y++) {
            for (long x = lattice.purge_ghost_count_x; x < lattice._size_sub_box_x + lattice.purge_ghost_count_x; x++) {
                if (_atoms[z][y][x].x[0] < domain.meas_global_box_coord_region.x_low - delta ||
                    _atoms[z][y][x].x[1] < domain.meas_global_box_coord_region.y_low - delta ||
                    _atoms[z][y][x].x[2] < domain.meas_global_box_coord_region.z_low - delta ||
                    _atoms[z][y][x].x[0] > domain.meas_global_box_coord_region.x_high + delta ||
                    _atoms[z][y][x].x[1] > domain.meas_global_box_coord_region.y_high + delta ||
                    _atoms[z][y][x].x[2] > domain.meas_global_box_coord_region.z_high + delta) {
                    return true;
                }
            }
        }
    }
    return false;
}

/**
 *******************************************
 * following (to end of this file) are the
 * implementations about atom list iterator.
 *******************************************
 */
AtomList::iterator AtomList::begin() {
    return AtomList::iterator(nullptr);
}

AtomList::iterator AtomList::end() {
    return AtomList::iterator(nullptr);
}

AtomList::const_iterator AtomList::begin() const {
    return AtomList::const_iterator(nullptr);
}

AtomList::const_iterator AtomList::end() const {
    return AtomList::const_iterator(nullptr);
}

// iterator
AtomList::iterator::iterator(AtomElement *ptr) : ptr_(ptr) {}

AtomList::iterator::self_type &AtomList::iterator::operator++() {
    ptr_++;
    return *this;
}

AtomList::iterator::self_type AtomList::iterator::operator++(int) {
    self_type i = *this;
    ptr_++;
    return i;
}

AtomElement &AtomList::iterator::operator*() {
    return *ptr_;
}

AtomElement *AtomList::iterator::operator->() {
    return ptr_;
}

bool AtomList::iterator::operator==(const AtomList::iterator::self_type &rhs) {
    return ptr_ == rhs.ptr_;
}

bool AtomList::iterator::operator!=(const AtomList::iterator::self_type &rhs) {
    return ptr_ != rhs.ptr_;
}

void AtomList::iterator::next() {

}

// const iterator
AtomList::const_iterator::const_iterator(AtomElement *ptr) : ptr_(ptr) {}

AtomList::const_iterator::self_type &AtomList::const_iterator::operator++() {
    ptr_++;
    return *this;
}

AtomList::const_iterator::self_type AtomList::const_iterator::operator++(int junk) {
    self_type i = *this;
    ptr_++;
    return i;
}

const AtomElement &AtomList::const_iterator::operator*() {
    return *ptr_;
}

bool AtomList::const_iterator::operator==(const AtomList::const_iterator::self_type &rhs) {
    return ptr_ == rhs.ptr_;
}

bool AtomList::const_iterator::operator!=(const AtomList::const_iterator::self_type &rhs) {
    return ptr_ != rhs.ptr_;
}
