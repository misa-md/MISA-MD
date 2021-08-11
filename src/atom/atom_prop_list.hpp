//
// Created by chu genshen on 2021/8/11.
//

#ifndef MISA_MD_ATOM_PROP_LIST_HPP
#define MISA_MD_ATOM_PROP_LIST_HPP

#include "arch/hardware_accelerate.hpp"

#include "../lattice/lattice.h"
#include "../types/pre_define.h"

template<typename T>
class AtomPropList {
private:
    T *data = nullptr;

    const BccLattice lattice;
public:
    explicit AtomPropList(const BccLattice &lattice) : lattice(lattice) {
        bool need_create = true;
        if (isArchAccSupport()) {
            // we may create memory using other api (e.g. pinned memory on CUDA platform).
            data = archCreateAtomsMemory(lattice._size_x, lattice._size_y, lattice._size_z);
            need_create = (data == nullptr);
        }
        if (need_create) {
            data = new AtomElement[lattice._size_z * lattice._size_x * lattice._size_y];
        }
    }

    inline T *_data() {
        return data;
    }

    void destroyPropList() {
        bool need_des = true;
        if (isArchAccSupport()) {
            need_des = !archReleaseAtomsMemory(data);
        }
        if (need_des) {
            delete[] data;
        }
    }


    /**
     * get atom element by linear index.
     * index = (zIndex * p_domain->getGhostLatticeSize(1) + yIndex) * p_domain->getGhostLatticeSize(0) + xIndex;
     * @param index
     * @return
     */
    inline T &getAtomEleByLinearIndex(_type_atom_index index) const {
        return data[index];
    }

    /**
     * calculate index of 1d atom array by the lattice coordinate in 3d
     * @param x, y, z the lattice coordinate in 3d
     * @return the index in 1d atoms array of the atom specified by @param x,y,z coordinate
     */
    inline _type_atom_index getAtomIndex(_type_atom_index x, _type_atom_index y, _type_atom_index z) const {
        return (z * lattice._size_y + y) * lattice._size_x + x;
    }

    /**
     * get the reference of {@class AtomElement} by the lattice index in sub-box(not include ghost lattice).
     * @param index_x lattice index in sub-box at x dimension.
     * @param index_y lattice index in sub-box at y dimension.
     * @param index_z lattice index in sub-box at z dimension.
     * @return reference of {@class AtomElement} at corresponding position.
     */
    inline T &
    getAtomEleBySubBoxIndex(_type_atom_index index_x, _type_atom_index index_y, _type_atom_index index_z) {
        return data[getAtomIndex(lattice.purge_ghost_count_x + index_x,
                                 lattice.purge_ghost_count_y + index_y,
                                 lattice.purge_ghost_count_z + index_z)];
    }

    /**
     * just by index of ghost lattice.
     * @param index_x
     * @param index_y
     * @param index_z
     * @return
     */
    inline T &
    getAtomEleByGhostIndex(_type_atom_index index_x, _type_atom_index index_y, _type_atom_index index_z) const {
        return data[getAtomIndex(index_x, index_y, index_z)];
    }
};

#endif //MISA_MD_ATOM_PROP_LIST_HPP
