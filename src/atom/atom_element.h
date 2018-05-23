//
// Created by genshen on 5/7/18.
//

#ifndef CRYSTAL_MD_ATOM_ELEMENT_H
#define CRYSTAL_MD_ATOM_ELEMENT_H


#include "../pre_define.h"
#include "atom_types.h"

/**
 * This class describes the attributes of one atom, such as location, force, velocity.etc.
 */

typedef atom_type::atom_type _type_atom_type_enum;

class AtomElement {
public:
    _type_atom_id id; // atom id.
    // <del> @deprecated
    _type_atom_type type; // atom type @depresed
    // </del>
    _type_atom_type_enum _tp; // atom type

    _type_atom_location x[DIMENSION]; // atom position.
    _type_atom_velocity v[DIMENSION]; // atom velocity.
    _type_atom_force f[DIMENSION];  // atom force.

    _type_atom_rho rho; // electron charge density
    _type_atom_df df; // embedded energy

    /**
     * check whether this atom is Inter atom.
     * @return
     */
    bool isInterElement() const;

};

#endif //CRYSTAL_MD_ATOM_ELEMENT_H
