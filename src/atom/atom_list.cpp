//
// Created by genshen on 5/8/18.
//

#include <utils/mpi_domain.h>
#include <pack/lat_particle_packer.h>
#include <comm.hpp>
#include "atom_list.h"
#include "../pack/pack.h"
#include "../utils/mpi_data_types.h"

AtomList::AtomList(_type_atom_count size_x, _type_atom_count size_y, _type_atom_count size_z,
                   _type_atom_count size_sub_box_x, _type_atom_count size_sub_box_y, _type_atom_count size_sub_box_z,
                   _type_atom_count ghost_count_x, _type_atom_count ghost_count_y, _type_atom_count ghost_count_z) :
        _size(size_x * size_y * size_z),
        _size_x(size_x), _size_y(size_y), _size_z(size_z),
        _size_sub_box_x(size_sub_box_x),
        _size_sub_box_y(size_sub_box_y),
        _size_sub_box_z(size_sub_box_z),
        purge_ghost_count_x(ghost_count_x),
        purge_ghost_count_y(ghost_count_y),
        purge_ghost_count_z(ghost_count_z) {

    _atoms = new AtomElement **[size_z];
    for (_type_atom_count z = 0; z < size_z; z++) {
        _atoms[z] = new AtomElement *[size_y];
        for (_type_atom_count y = 0; y < size_y; y++) {
            _atoms[z][y] = new AtomElement[size_x];
        }
    }
}

AtomList::~AtomList() {
    for (_type_atom_count z = 0; z < _size_z; z++) {
        for (_type_atom_count y = 0; y < _size_y; y++) {
            delete[] _atoms[z][y];
        }
        delete[] _atoms[z];
    }
    delete[] _atoms;
}

void AtomList::exchangeAtomFirst(comm::Domain *p_domain, int cutlattice) {
    sendlist.resize(6);
    recvlist.resize(6);
    for (unsigned short d = 0; d < DIMENSION; d++) {
        for (int direction = LOWER; direction <= HIGHER; direction++) {
            // 找到要发送给邻居的原子
            switch (d) {
                case 0:
                    getatomx(p_domain, cutlattice, direction, sendlist);
                    break;
                case 1:
                    getatomy(p_domain, cutlattice, direction, sendlist);
                    break;
                case 2:
                    getatomz(p_domain, cutlattice, direction, sendlist);
                    break;
                default:
                    break;
            }
        }
    }

    LatPackerFirst lat_packer(*p_domain, *this, sendlist, recvlist);
    comm::neiSendReceive<LatParticleData>(&lat_packer,
                                          MPIDomain::toCommProcess(),
                                          mpi_types::_mpi_latParticle_data,
                                          p_domain->rank_id_neighbours);
}

void AtomList::exchangeAtom(comm::Domain *p_domain) {
    LatPacker lat_packer(*p_domain, *this, sendlist, recvlist);
    comm::neiSendReceive<LatParticleData>(&lat_packer,
                                          MPIDomain::toCommProcess(),
                                          mpi_types::_mpi_latParticle_data,
                                          p_domain->rank_id_neighbours);
}

void AtomList::getatomx(comm::Domain *p_domain, int _cutlattice, int direction,
                        std::vector<std::vector<_type_atom_id> > &sendlist) {
    _type_atom_id i;
    if (direction == 0) {
        //找到要发送到邻居进程的区域
        int xstart = p_domain->dbx_lattice_size_ghost[0];
        int ystart = p_domain->dbx_lattice_size_ghost[1];
        int zstart = p_domain->dbx_lattice_size_ghost[2];
        int xstop = xstart + p_domain->dbx_lattice_size_ghost[0]; // note: this is ghost lattice size.
        int ystop = ystart + p_domain->dbx_lattice_size_sub_box[1];
        int zstop = zstart + p_domain->dbx_lattice_size_sub_box[2];

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[0].push_back(i);
                }
            }
        }
    } else {
        //找到要发送到邻居进程的区域
        int xstart = p_domain->dbx_lattice_size_ghost[0] + p_domain->dbx_lattice_size_sub_box[0] - ((_cutlattice) * 2);
        int ystart = p_domain->dbx_lattice_size_ghost[1];
        int zstart = p_domain->dbx_lattice_size_ghost[2];
        int xstop = p_domain->dbx_lattice_size_ghost[0] + p_domain->dbx_lattice_size_sub_box[0];
        int ystop = ystart + p_domain->dbx_lattice_size_sub_box[1];
        int zstop = zstart + p_domain->dbx_lattice_size_sub_box[2];

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[1].push_back(i);
                }
            }
        }
    }
}

void AtomList::getatomy(comm::Domain *p_domain, int _cutlattice, int direction,
                        std::vector<std::vector<_type_atom_id> > &sendlist) {
    int i;
    if (direction == 0) {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = p_domain->dbx_lattice_size_ghost[1];
        int zstart = p_domain->dbx_lattice_size_ghost[2];
        int xstop = p_domain->dbx_lattice_size_ghost_extended[0];
        int ystop = ystart + _cutlattice;
        int zstop = zstart + p_domain->dbx_lattice_size_sub_box[2];

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[2].push_back(i);
                }
            }
        }
    } else {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = p_domain->dbx_lattice_size_ghost[1] + p_domain->dbx_lattice_size_sub_box[1] - (_cutlattice);
        int zstart = p_domain->dbx_lattice_size_ghost[2];
        int xstop = p_domain->dbx_lattice_size_ghost_extended[0];
        int ystop = p_domain->dbx_lattice_size_ghost[1] + p_domain->dbx_lattice_size_sub_box[1];
        int zstop = zstart + p_domain->dbx_lattice_size_sub_box[2];

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[3].push_back(i);
                }
            }
        }
    }
}

void AtomList::getatomz(comm::Domain *p_domain, int _cutlattice, int direction,
                        std::vector<std::vector<_type_atom_id> > &sendlist) {
    int i;
    if (direction == 0) {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = 0;
        int zstart = p_domain->dbx_lattice_size_ghost[2];
        int xstop = p_domain->dbx_lattice_size_ghost_extended[0];
        int ystop = p_domain->dbx_lattice_size_ghost_extended[1];
        int zstop = zstart + _cutlattice;

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[4].push_back(i);
                }
            }
        }
    } else {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = 0;
        int zstart = p_domain->dbx_lattice_size_ghost[2] + p_domain->dbx_lattice_size_sub_box[2] - (_cutlattice);
        int xstop = p_domain->dbx_lattice_size_ghost_extended[0];
        int ystop = p_domain->dbx_lattice_size_ghost_extended[1];
        int zstop = p_domain->dbx_lattice_size_ghost[2] + p_domain->dbx_lattice_size_sub_box[2];

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[5].push_back(i);
                }
            }
        }
    }
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
