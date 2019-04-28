//
// Created by genshen on 2019-04-11.
//

template<class T>
NeighbourIndex<T>::NeighbourIndex(AtomList &atom_list)
        :atom_list(atom_list), nei_even_offsets(), nei_odd_offsets(),
         nei_half_even_offsets(), nei_half_odd_offsets() {}

template<class T>
void NeighbourIndex<T>::make(const _type_lattice_size cut_lattice,
                             const double cutoff_radius_factor) {
    const double cutoff_lat_factor = cutoff_radius_factor + 1.0 / 2 + 1.0 / 2;
    // if x index of a particle is even (the particle is lattice point,晶格点).
    for (_type_atom_index zIndex = -cut_lattice;
         zIndex <= cut_lattice; zIndex++) { // loop for (2*_cutlattice + 1) times.
        for (_type_atom_index yIndex = -cut_lattice; yIndex <= cut_lattice; yIndex++) {
            for (_type_atom_index xIndex = -cut_lattice; xIndex <= cut_lattice; xIndex++) {
//               uint z = (double) zIndex + (((double) (xIndex % 2)) / 2); // zIndex plus 1/2 (odd) or 0(even).
//               uint y = (double) yIndex + (((double) (xIndex % 2)) / 2);
//               uint x = (double) xIndex / 2;
                // lattice neighbour points whose index is (2*Index,yIndex,zIndex).
                {
                    const double r = xIndex * xIndex + yIndex * yIndex + zIndex * zIndex;
                    // r > 0 means neighbour index can not be itself.
                    if (r < cutoff_lat_factor * cutoff_lat_factor && r > 0) {
                        nei_even_offsets.push_back(NeiOffset{2 * xIndex, yIndex, zIndex});
                        if (isPositiveIndex(xIndex, yIndex, zIndex)) {
                            nei_half_even_offsets.push_back(NeiOffset{2 * xIndex, yIndex, zIndex});
                        }
                    }
                }
                // bcc body center neighbour points whose index is (2*Index+1,yIndex,zIndex).
                {
                    const double r = (xIndex + 0.5) * (xIndex + 0.5) +
                                     (yIndex + 0.5) * (yIndex + 0.5) +
                                     (zIndex + 0.5) * (zIndex + 0.5);
                    if (r < cutoff_lat_factor * cutoff_lat_factor && r > 0) { // in fact "r > 0" is not used.
                        nei_even_offsets.push_back(NeiOffset{2 * xIndex + 1, yIndex, zIndex});
                        if (isPositiveIndex(xIndex + 0.5, yIndex + 0.5, zIndex + 0.5)) {
                            nei_half_even_offsets.push_back(NeiOffset{2 * xIndex + 1, yIndex, zIndex});
                        }
                    }
                }
            }
        }
    }
    // if x index of a particle is odd (the particle is BCC body center point,体心).
    for (_type_atom_index zIndex = -cut_lattice;
         zIndex <= cut_lattice; zIndex++) { // loop for (2*_cutlattice + 1) times.
        for (_type_atom_index yIndex = -cut_lattice; yIndex <= cut_lattice; yIndex++) {
            for (_type_atom_index xIndex = -cut_lattice; xIndex <= cut_lattice; xIndex++) {
                // BCC body center neighbour points whose index is (2*Index,yIndex,zIndex).
                {
                    const double r = xIndex * xIndex + yIndex * yIndex + zIndex * zIndex;
                    if (r < cutoff_lat_factor * cutoff_lat_factor && r > 0) {
                        nei_odd_offsets.push_back(NeiOffset{2 * xIndex, yIndex, zIndex});
                        if (isPositiveIndex(xIndex, yIndex, zIndex)) {
                            nei_half_odd_offsets.push_back(NeiOffset{2 * xIndex, yIndex, zIndex});
                        }
                    }
                }
                // lattice neighbour points whose index is (2*Index-1,yIndex,zIndex).
                {
                    const double r = (xIndex - 0.5) * (xIndex - 0.5) +
                                     (yIndex - 0.5) * (yIndex - 0.5) +
                                     (zIndex - 0.5) * (zIndex - 0.5);
                    if (r < cutoff_lat_factor * cutoff_lat_factor && r > 0) {
                        nei_odd_offsets.push_back(NeiOffset{2 * xIndex - 1, yIndex, zIndex});
                        if (isPositiveIndex(xIndex - 0.5, yIndex - 0.5, zIndex - 0.6)) {
                            nei_half_odd_offsets.push_back(NeiOffset{2 * xIndex - 1, yIndex, zIndex});
                        }
                    }
                }
            }
        }
    }
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
    const NeiOffset offset{x, y, z};
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
    const NeiOffset offset{x, y, z};
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
