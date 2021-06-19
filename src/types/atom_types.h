//
// Created by genshen on 2018-05-06.
//

#ifndef MISA_MD_ATOM_TYPES_H
#define MISA_MD_ATOM_TYPES_H

#include <cstdio>
#include <vector>
#include <string>

/**
 * @see also config_values.h#AlloyRatio
 */
#include "pre_define.h"

namespace atom_type {

    enum atom_type { // atom id
        INVALID = -1 // Fe = 0, Cu = 1, Ni = 2 /*,Co */
    };

    /**
     * mass array of system atoms.
     * atom id is the index of this array.
     */
    extern std::vector<double> mass_array;
    /**
     * the count of atom types.
     * How many types of atoms appears in this system.
     * private variable.
     */
    extern _type_atom_types num_atom_types;

    void setGlobalAtomMasses(const std::vector<double> masses);

    void releaseGlobalAtomMesses();

    /**
     * get element relative atomic mass.
     * @param type
     * @return relative atomic mass
     */
    inline _type_atom_mass getAtomMass(atom_type tp) {
        if (tp < 0 || tp >= num_atom_types) {
            return 0.0;
        }
        return mass_array[tp];
    }

    /**
     * get the ith atom type.
     * @param i index of atom numbered starting 0.
     * @return atom_type.
     */
    inline atom_type getAtomTypeByNum(int i) {
        return static_cast<atom_type>(i);
    }

    // tod return type atom_type::_type_prop_key
    inline unsigned short getTypeIdByType(atom_type tp) {
        if (tp < 0 || tp >= num_atom_types) {
            printf("waring, not expect id zero, it may cause error.\n");
            return 0;
        }
        switch (tp) {
            case 0:
                return 0;
            case 1:
                return 1;
            case 2:
                return 2;
            default:
                return 0;
        }
        return tp;
    }

//    inline int getAtomTypeIndex(atom_type type) {
//        return type;
//    }
//
//    inline int getAtomTypeMaterialsIndex(atom_type type1, atom_type type2) {
//        return num_atom_types * type1 + type2;
//    }
}

#endif //MISA_MD_ATOM_TYPES_H
