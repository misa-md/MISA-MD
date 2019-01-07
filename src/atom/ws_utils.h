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

    /**
     * It returns reference of a atom who has a least distance
     * from its lattice coordinate to the atom(source atom) specified by @param atom.
     * @note make sure the source atom is in the sub-box or in ghost area before
     * calling this method (this method does not guarantee that).
     * @param atom_list the atom list.
     * @param atom the source atom.
     * @param p_domain box domain.
     * @return an atom reference who has a least distance from its lattice coordinate to the source atom.
     */
    AtomElement &findNearLatAtom(AtomList *atom_list, const AtomElement &src_atom, const Domain *p_domain);

    /**
     * find an near atom from all lattice atoms in sub box (the ghost area in not included).
     * The near atom is defined as: the atom has a least distance from
     * its lattice coordinate to the atom(source atom) specified by @param atom.
     *
     * If there is no near atom in sub box found, an nullptr pointer will be returned.
     * @param atom_list pointer to the atom list.
     * @param src_atom reference to source atom.
     * @param p_domain pointer to box domain.
     * @return pointer to the near atom or nullptr if the near atom is not found.
     */
    AtomElement *finNearLatAtomInSubBox(AtomList *atom_list, const AtomElement &src_atom, const Domain *p_domain);
};


#endif //CRYSTALMD_WS_UTILS_H
