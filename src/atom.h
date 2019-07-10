//
// Created by baihe back to 2016-01-06.
//

#ifndef CRYSTAL_MD_ATOM_H
#define CRYSTAL_MD_ATOM_H

#include <cstdio>
#include <vector>
#include <eam.h>

#include <domain/bcc_domain.h>
#include <logs/logs.h>

#include "atom/atom_set.h"
#include "atom/atom_element.h"
#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"
#include "pack/particledata.h"
#include "pack/lat_particle_data.h"

class atom : public AtomSet {
public :
    atom(comm::BccDomain *domain);

    /**
     * move atoms to inter-atom list if the atoms is not in its lattice.
     * @return n_flag
     */
    int decide();

    void clearForce();

    void computeEam(eam *pot, double &comm);

    /**
     * set velocity of a atom whose position is specified by array @param lat
     * This atom is called PKA (Primary Knock-on Atom).
     * We first compute the the velocity in each direction from the energy, and the set velocity of the atom.
     * The unit of energy is ev, the unit of velocity is 100m/s.
     *
     * @param lat the position of pka.
     * @param direction the vector of velocity to be set.
     * @param energy the energy of PKA.
     */
    void setv(int lat[4], double direction[3], double energy);


    void sendForce();

private:
    comm::BccDomain *p_domain;

    /**
     * calculate electron density for all lattice atoms.
     * @param pot pointer of eam potential object.
     */
    void latRho(eam *pot);

    /**
     * calculate electron density for all interstitial atoms.
     * @param pot pointer of eam potential object.
     */
    void interRho(eam *pot);

    /**
     * calculate derivative of embedded energy for all lattice atoms.
     * @param pot pointer of eam potential object.
     */
    void latDf(eam *pot);

    /**
     * calculate force for all lattice atoms.
     * @param pot pointer of eam potential object.
     */
    void latForce(eam *pot);

    /**
     * calculate force for all interstitial atoms.
     * @param pot pointer of eam potential object.
     */
    void interForce(eam *pot);
};

#endif // CRYSTAL_MD_ATOM_H
