//
// Created by genshen on 5/7/18.
//

#ifndef MISA_MD_ATOM_ELEMENT_H
#define MISA_MD_ATOM_ELEMENT_H


#include "../types/pre_define.h"
#include "../types/atom_types.h"

/**
 * This class describes the attributes of one atom, such as location, force, velocity.etc.
 */

typedef atom_type::atom_type _type_atom_type_enum;

class AtomElement {
public:
    _type_atom_id id; // atom id.
    // <del> @deprecated
//    _type_atom_type type; // atom type @depressed
    // </del>
    _type_atom_type_enum type; // atom type

    _type_atom_location x[DIMENSION]; // atom position.
    _type_atom_velocity v[DIMENSION]; // atom velocity.
    _type_atom_force f[DIMENSION];  // atom force.

    _type_atom_rho rho; // electron charge density
    _type_atom_df df; // embedded energy

    /**
     * check whether this atom is Inter atom.
     * @return
     */
    inline bool isInterElement() const {
        return type == atom_type::INVALID;
    }

};

#endif //MISA_MD_ATOM_ELEMENT_H
