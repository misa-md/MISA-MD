//
// Created by genshen on 2019-01-02.
//

#ifndef MISA_MD_WS_UTILS_H
#define MISA_MD_WS_UTILS_H

#include <comm/domain/domain.h>
#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"
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
     * @param src_atom position of the source atom
     * @param p_domain the domain of current processor
     * @return the flag of status.
     */
    const box::_type_flag_32 isOutBox(_type_atom_location src_atom_x[DIMENSION], const comm::Domain *p_domain);

    /**
     * It returns an array index who has a least distance
     * from its lattice coordinate to the atom(source atom) specified by @param atom.
     * @note make sure the source atom is in the sub-box or in ghost area before
     * calling this method (this method does not guarantee that).
     * @deprecated use findNearLatIndexInSubBox instead.
     *
     * @param atom_list the atom list.
     * @param atom the source atom.
     * @param p_domain box domain.
     * @return an array index who has a least distance from its lattice coordinate to the source atom.
     */
    _type_atom_index findNearLatAtom(AtomList *atom_list, const AtomElement &src_atom, const comm::Domain *p_domain);

    /**
     * find an near atom from all lattice atoms in sub box (the ghost area in not included).
     * The nearest atom is defined as: the atom has a least distance from
     * its lattice coordinate to the atom(source atom) specified by @param atom.
     * @deprecated use findNearLatIndexInSubBox instead.
     *
     * If there is no nearest atom in sub box found, an nullptr pointer will be returned.
     * @param atom_list pointer to the atom list.
     * @param src_atom reference to source atom.
     * @param p_domain pointer to box domain.
     * @return array index to the nearest atom or nullptr if the nearest atom is not found.
     */
    _type_atom_index
    findNearLatAtomInSubBox(AtomList *atom_list, const AtomElement &src_atom, const comm::Domain *p_domain);

    /**
     * It returns the lattice index in 3d(not atom) of nearest lattice of @param src_atom.
     *
     * if the nearest atom is not found in sub-box(the @param src_atom can be out of sub-box),
     * box::IndexNotExists will be returned.
     *
     * @param lattice the lattice which described the simulation box.
     * @param src_atom reference to source atom.
     * @param p_domain pointer to box domain.
     * @return lattice index of nearest atom in 3d, or box::IndexNotExists if the nearest atom is not found
     */
    _type_atom_index findNearLatIndexInSubBox(const BccLattice &lattice, const AtomElement &src_atom, const comm::Domain *p_domain);

    /**
     * get the coordinate(starting from ghost area, not sub-box) of nearest lattice of @param src_atom.
     *
     * @note the x lattice coordinate is doubled.
     * @param src_atom source atom.
     * @param p_domain pointer to box domain.
     * @param coords the coordinate of nearest lattice to be returned.
     */
    void getNearLatCoord(const AtomElement &src_atom, const comm::Domain *p_domain, _type_atom_index coords[DIMENSION]);

    /**
     * similar as above, but the lattice coordinate is relative to sub-box, ghost area is not included.
     *
     * @note the x lattice coordinate is doubled.
     * @param src_atom source atom.
     * @param p_domain pointer to box domain.
     * @param coords the coordinate of nearest lattice to be returned.
     */
    void getNearLatSubBoxCoord(const AtomElement &src_atom, const comm::Domain *p_domain, _type_atom_index coords[DIMENSION]);

    /**
     * check weather an atom's coordinate is in current sub-box (the sub box on current processor).
     *
     * In the implementations, we first calculate the nearest lattice coordinate of @param src_atom
     * and check whether the lattice coordinate is in the sub-box.
     * We does not check the real coordinate of source atom (it is not very elegant).
     * @return true for in the sub-box, false for otherwise.
     */
    bool isInBox(const AtomElement &src_atom, const comm::Domain *p_domain);
};


#endif //MISA_MD_WS_UTILS_H
