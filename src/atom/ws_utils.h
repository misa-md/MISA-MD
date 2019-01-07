//
// Created by genshen on 2019-01-02.
//

#ifndef CRYSTALMD_WS_UTILS_H
#define CRYSTALMD_WS_UTILS_H

#include "atom_list.h"
#include "inter_atom_list.h"
#include "box.h"

/**
 * utils for Wigner Seitz  cell.
 *
 * The 3D-view of the Wigner-Seitz cell for bcc-structure can be found at:
 * https://www.physics-in-a-nutshell.com/article/12/body-centered-cubic-bcc
 */
namespace ws {
    /**
     * check weather an atom is in the sub-box of current processor.
     * @param src_atom the referent of atom
     * @param p_domain the domain of current processor
     * @return the flag of status.
     */
    const box::_type_flag_32 isOutBox(const AtomElement &src_atom, const Domain *p_domain);
};


#endif //CRYSTALMD_WS_UTILS_H
