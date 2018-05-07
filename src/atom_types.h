//
// Created by genshen on 2018-05-06.
//

#ifndef CRYSTAL_MD_ATOM_TYPES_H
#define CRYSTAL_MD_ATOM_TYPES_H

/**
 * @see also config_values.h#AlloyRatio
 */

#define _ATOM_TYPES 3 // 3 different types of atoms(Fe-Cu-Ni).

// WARNING: DO NOT use the marco below directly.
#define _R_A_M_Fe 55.845 // in which, R_A_M means "relative atomic mass"
#define _R_A_M_Cu 63.546
#define _R_A_M_Ni 58.6934

namespace atom_type {
    typedef short _type_atom_types;

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
    inline double getAtomMass(atom_type atom) {
        switch (atom) {
            case Fe:
                return _R_A_M_Fe;
            case Cu:
                return _R_A_M_Cu;
            case Ni:
                return _R_A_M_Ni;
        }
    }

}

#endif //CRYSTAL_MD_ATOM_TYPES_H
