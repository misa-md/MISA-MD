//
// Created by baihe back to 2016-01-06.
//

#ifndef CRYSTAL_MD_ATOM_H
#define CRYSTAL_MD_ATOM_H

#include <cstdio>
#include <vector>
#include <eam.h>

#include <domain/domain.h>

#include "atom/atom_set.h"
#include "atom/atom_element.h"
#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"
#include "pack/particledata.h"
#include "pack/lat_particle_data.h"

class atom : public AtomSet {
public :
    atom(comm::Domain *domain);

    /**
     * move atoms to inter-atom list if the atoms is not in its lattice.
     * @return n_flag
     */
    int decide();

    void clearForce();

    void computeEam(eam *pot, double &comm);

    void latRho(eam *pot, double &comm);

    void interRho(eam *pot, double &comm);

    void latDf(eam *pot, double &comm);

    void latForce(eam *pot, double &comm);

    void interForce(eam *pot, double &comm);

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

    void print_force(const std::string filename);

    void sendForce();

private:
    comm::Domain *p_domain;

};

#endif // CRYSTAL_MD_ATOM_H
