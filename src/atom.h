//
// Created by baihe back to 2016-01-06.
//

#ifndef MISA_MD_ATOM_H
#define MISA_MD_ATOM_H

#include <cstdio>
#include <vector>
#include <eam.h>

#include <comm/domain/bcc_domain.h>
#include <logs/logs.h>

#include "atom/atom_set.h"
#include "atom/atom_element.h"
#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"
#include "pack/send_recv_lists.h"
#include "pack/particledata.h"
#include "pack/lat_particle_data.h"

class atom : public AtomSet {
public :
    atom(comm::BccDomain *domain);

    ~ atom();

    /**
     * move atoms to inter-atom list if the atoms is not in its lattice.
     * @return n_flag
     */
    int decide();

    void clearForce();

    void computeEam(const unsigned short potentialType,eam *pot, double &comm);

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
    void setv(const _type_lattice_coord lat[4], const double direction[3], const double energy);

public:
    SendRecvLists *p_send_recv_list;

    /**
     * set velocity for a atom whose position is specified by array @param lat.
     * @param lat_x, lat_y, lat_z lattice position of the atom (x dimension is not doubled), whose velocity will be changed.
     * @param v the velocity value at x, y, z dimension.
     */
    void setv(const _type_lattice_coord lat_x, const _type_lattice_coord lat_y,
              const _type_lattice_coord lat_z, const double v[DIMENSION]);


private:
    comm::BccDomain *p_domain;

    /**
     * calculate electron density for all lattice atoms.
     * @param potentialType type of potentail calculation.
     * @param pot pointer of eam potential object.
     */
    void latRho(const unsigned short potentialType,eam *pot);

    /**
     * calculate electron density for all interstitial atoms.
     * @param potentialType type of potentail calculation.
     * @param pot pointer of eam potential object.
     */
    void interRho(const unsigned short potentialType,eam *pot);

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

#endif // MISA_MD_ATOM_H
