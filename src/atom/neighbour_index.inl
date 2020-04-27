//
// Created by genshen on 2019-04-11.
//
#include <cassert>

template<class T>
NeighbourIndex<T>::NeighbourIndex(AtomList &atom_list)
        :atom_list(atom_list), nei_even_offsets(), nei_odd_offsets(),
         nei_half_even_offsets(), nei_half_odd_offsets() {}

template<class T>
void NeighbourIndex<T>::make(const _type_lattice_size cut_lattice,
                             const double cutoff_radius_factor) {
    const double cutoff_lat_factor = cutoff_radius_factor + 1.0 / 2 + 1.0 / 2;
    // if x index of a particle is even (the particle is lattice point,晶格点).
    for (_type_atom_index zIndex = -cut_lattice - 1;
         zIndex <= cut_lattice + 1; zIndex++) { // loop for (2*_cutlattice + 1) times.
        for (_type_atom_index yIndex = -cut_lattice - 1; yIndex <= cut_lattice + 1; yIndex++) {
            for (_type_atom_index xIndex = -2 * cut_lattice - 2; xIndex <= 2 * cut_lattice + 2; xIndex++) {
                // lattice neighbour points whose index is (2*Index,yIndex,zIndex).
                double z = (double) zIndex + (((double) (xIndex % 2)) / 2); // zIndex plus 1/2 (odd) or 0(even).
                double y = (double) yIndex + (((double) (xIndex % 2)) / 2);
                double x = ((double) xIndex) / 2;
                const double r = x * x + y * y + z * z;

                // r > 0 means neighbour index can not be itself.
                if (r < cutoff_lat_factor * cutoff_lat_factor && r > 0) {
                    // if xIndex is less than 0 and it is odd.
                    // for example xIndex=-1,yIndex=0,zIndex=0 =>
                    // then (x,y,z) is (-0.5,-0.5,-0.5) and (ix,iy,iz) = (-1,-1,-1) not (-1,0,0)
                    const _type_atom_index ix = xIndex;
                    const _type_atom_index iy = (xIndex < 0 && xIndex % 2 != 0) ? yIndex - 1 : yIndex;
                    const _type_atom_index iz = (xIndex < 0 && xIndex % 2 != 0) ? zIndex - 1 : zIndex;
                    const _type_atom_index even_offset = atom_list.lattice.IndexOf3DIndex(ix, iy, iz);
                    nei_even_offsets.push_back(even_offset);
                    if (isPositiveIndex(x, y, z)) {
                        nei_half_even_offsets.push_back(even_offset);
                    }
                }
            }
        }
    }
    // if x index of a particle is odd (the particle is BCC body center point,体心).
    for (_type_atom_index zIndex = -cut_lattice - 1;
         zIndex <= cut_lattice + 1; zIndex++) { // loop for (2*_cutlattice + 1) times.
        for (_type_atom_index yIndex = -cut_lattice - 1; yIndex <= cut_lattice + 1; yIndex++) {
            for (_type_atom_index xIndex = -2 * cut_lattice - 2; xIndex <= 2 * cut_lattice + 2; xIndex++) {
                // BCC body center neighbour points whose index is (2*Index,yIndex,zIndex).
                double z = (double) zIndex - (((double) (xIndex % 2)) / 2);
                double y = (double) yIndex - (((double) (xIndex % 2)) / 2);
                double x = (double) xIndex / 2;
                const double r = x * x + y * y + z * z;

                if (r < cutoff_lat_factor * cutoff_lat_factor && r > 0) {
                    // if xIndex,yIndex,zIndex is (-1,0,0), then (x,y,z) is (-0.5,0.5,0.5)
                    // and (ix,iy,iz) is (-1,1,1)
                    const _type_atom_index ix = xIndex;
                    const _type_atom_index iy = (xIndex < 0 && xIndex % 2 != 0) ? yIndex + 1 : yIndex;
                    const _type_atom_index iz = (xIndex < 0 && xIndex % 2 != 0) ? zIndex + 1 : zIndex;
                    const _type_atom_index odd_offset = atom_list.lattice.IndexOf3DIndex(ix, iy, iz);
                    nei_odd_offsets.push_back(odd_offset);
                    if (isPositiveIndex(x, y, z)) {
                        nei_half_odd_offsets.push_back(odd_offset);
                    }
                }
            }
        }
    }
#ifdef MD_DEV_MODE
    assert(nei_odd_offsets.size() == 2 * nei_half_odd_offsets.size());
    assert(nei_even_offsets.size() == 2 * nei_half_even_offsets.size());
    assert(nei_odd_offsets.size() == nei_even_offsets.size());
#endif
}

template<class T>
bool NeighbourIndex<T>::isPositiveIndex(const double x, const double y, const double z) {
    if (z > 0) { // todo float equal
        return true;
    } else if (z == 0) {
        if (y > 0) {
            return true;
        } else if (y == 0) {
            if (x > 0) {
                return true;
            }
        }
    }
    return false;
}

template<class T>
NeiIterator<T, T &, T *>
NeighbourIndex<T>::begin(const bool half_itl, const _type_atom_index x,
                         const _type_atom_index y, const _type_atom_index z) {
    const NeiOffset offset = atom_list.lattice.IndexOf3DIndex(x, y, z);
    const int flag = (half_itl ? 2 : 0) | (x % 2 == 0 ? 1 : 0);
    switch (flag) {
        case 0: // 0b00
            return NeiIterator<T, T &, T *>(&nei_odd_offsets, offset, &atom_list);
        case 1: //  0b01
            return NeiIterator<T, T &, T *>(&nei_even_offsets, offset, &atom_list);
        case 2: // 0b10
            return NeiIterator<T, T &, T *>(&nei_half_odd_offsets, offset, &atom_list);
        case 3: // 0b11
            return NeiIterator<T, T &, T *>(&nei_half_even_offsets, offset, &atom_list);
        default: // default is not used.
            return NeiIterator<T, T &, T *>(&nei_even_offsets, offset, &atom_list);
    }
}

template<class T>
NeiIterator<T, T &, T *>
NeighbourIndex<T>::end(const bool half_itl, const _type_atom_index x,
                       const _type_atom_index y, const _type_atom_index z) {
    const NeiOffset offset = atom_list.lattice.IndexOf3DIndex(x, y, z);
    const int flag = (half_itl ? 2 : 0) | (x % 2 == 0 ? 1 : 0);
    switch (flag) {
        case 0: // 0b00
            return NeiIterator<T, T &, T *>(&nei_odd_offsets, offset, &atom_list, nei_odd_offsets.size());
        case 1: //  0b01
            return NeiIterator<T, T &, T *>(&nei_even_offsets, offset, &atom_list, nei_even_offsets.size());
        case 2: // 0b10
            return NeiIterator<T, T &, T *>(&nei_half_odd_offsets, offset, &atom_list, nei_half_odd_offsets.size());
        case 3: // 0b11
            return NeiIterator<T, T &, T *>(&nei_half_even_offsets, offset, &atom_list, nei_half_even_offsets.size());
        default: // default is not used.
            return NeiIterator<T, T &, T *>(&nei_even_offsets, offset, &atom_list, nei_even_offsets.size());
    }
}
