//
// Created by genshen on 2018-05-06.
//

#ifndef CRYSTAL_MD_ATOM_TYPES_H
#define CRYSTAL_MD_ATOM_TYPES_H

/**
 * @see also config_values.h#AlloyRatio
 */
#include "../pre_define.h"

#define _ATOM_TYPES 3 // 3 different types of atoms(Fe-Cu-Ni).

// WARNING: DO NOT use the marco below directly.
#define _R_A_M_Fe 55.845 // in which, R_A_M means "relative atomic mass"
#define _R_A_M_Cu 63.546
#define _R_A_M_Ni 58.6934

namespace atom_type {

    enum atom_type {
        Fe, Cu, Ni /*,Co */
    };

    /**
     * the count of atom types.
     * How many types of atoms appears in this system.
     */
    const _type_atom_types num_atom_types = _ATOM_TYPES;

    /**
     * get element relative atomic mass.
     * @param type
     * @return relative atomic mass
     */
    inline _type_atom_mass getAtomMass(atom_type atom) {
        switch (atom) {
            case Fe:
                return _R_A_M_Fe;
            case Cu:
                return _R_A_M_Cu;
            case Ni:
                return _R_A_M_Ni;
        }
    }

    /**
     * get the ith atom type.
     * @param i index of atom numbered starting 0.
     * @return atom_type.
     */
    inline atom_type getAtomTypeByNum(int i) {
        return static_cast<atom_type>(i);
    }

//    inline int getAtomTypeIndex(atom_type type) {
//        return type;
//    }
//
//    inline int getAtomTypeMatrisIndex(atom_type type1, atom_type type2) {
//        return num_atom_types * type1 + type2;
//    }
}

#endif //CRYSTAL_MD_ATOM_TYPES_H