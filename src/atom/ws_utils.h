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
     * It returns reference of a atom(nearest atom) who has a least distance
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
     * The nearest atom is defined as: the atom has a least distance from
     * its lattice coordinate to the atom(source atom) specified by @param atom.
     *
     * If there is no nearest atom in sub box found, an nullptr pointer will be returned.
     * @param atom_list pointer to the atom list.
     * @param src_atom reference to source atom.
     * @param p_domain pointer to box domain.
     * @return pointer to the nearest atom or nullptr if the nearest atom is not found.
     */
    AtomElement *findNearLatAtomInSubBox(AtomList *atom_list, const AtomElement &src_atom, const Domain *p_domain);

    /**
     * It returns the lattice index in 3d(not atom) of nearest lattice of @param src_atom.
     *
     * if the nearest atom is not found in sub-box(the @param src_atom can be out of sub-box),
     * box::IndexNotExists will be returned.
     *
     * @param atom_list pointer to the atom list.
     * @param src_atom reference to source atom.
     * @param p_domain pointer to box domain.
     * @return lattice index of nearest atom in 3d, or box::IndexNotExists if the nearest atom is not found
     */
    _type_atom_index findNearLatIndexInSubBox(AtomList *atom_list, const AtomElement &src_atom, const Domain *p_domain);

    /**
     * get the coordinate(starting from ghost area, not sub-box) of nearest lattice of @param src_atom.
     *
     * @note the x lattice coordinate is doubled.
     * @param src_atom source atom.
     * @param p_domain pointer to box domain.
     * @param coords the coordinate of nearest lattice to be returned.
     */
    void getNearLatCoord(const AtomElement &src_atom, const Domain *p_domain, _type_atom_index coords[DIMENSION]);

    /**
     * similar as above, but the lattice coordinate is relative to sub-box, ghost area is not included.
     *
     * @note the x lattice coordinate is doubled.
     * @param src_atom source atom.
     * @param p_domain pointer to box domain.
     * @param coords the coordinate of nearest lattice to be returned.
     */
    void getNearLatSubBoxCoord(const AtomElement &src_atom, const Domain *p_domain, _type_atom_index coords[DIMENSION]);
};


#endif //CRYSTALMD_WS_UTILS_H
