//
// Created by baihe back to 2016-12-22.
//

#ifndef CRYSTAL_MD_DOMAIN_DECOMPOSITION_H
#define CRYSTAL_MD_DOMAIN_DECOMPOSITION_H

#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR

#include <mpi.h>
#include <vector>

#include "atom.h"
#include "pre_config.h"

#define LOWER  0
#define HIGHER 1

class atom;

/**
 * If N can be decomposed as N = N_x * N_y * N_z, where N, N_x, N_y, N_z are all integer bigger than or equal to 1,
 * then the whole simulation box will be divided into N sub-box with N_z levels in z axis,
 * and in each level, it has  N_x * N_y sub-boxes.
 * Then, we can bind each processor to a cartesian coordinate (x,y,z) due to the boxes partition,
 * where 0 <= x < N_x, 0 <= z < N_z, 0 <= z < N_z.
 * Last, based on the cartesian coordinate (x,y,z),
 * each processor can get the cartesian coordinate of its contiguous sub-boxes.
 */
class DomainDecomposition {
public:
    friend class atom;

    /**
     * global information for the simulation box.
     */
    double _globalLength[DIMENSION]; // todo private.

    /**
     * the lower and upper bound coordinate for the global simulation box.
     */
    double _coord_global_box_low[DIMENSION], _coord_global_box_high[DIMENSION];

    DomainDecomposition();

    ~DomainDecomposition();

    /**
     * get global length of the simulation box at dimension {@var d}.
     * @param d dimension
     * @return box length at dimension d.
     */
    inline double getGlobalLength(int d) const {
        return _globalLength[d];
    }

    /**
     * set length of global simulation box.
     * and set upper and lower bound of global simulation box.
     * @param phaseSpace  phase space, the unit length of simulation box.
     * @param latticeConst lattice const
     */
    void establishGlobalDomain(const int64_t phaseSpace[DIMENSION], const double latticeConst);

    /**
     * set bound for current sub-box.
     * @param phaseSpace  phase space, the unit length of simulation box.
     * @param latticeConst lattice const
     */
    void establishLocalBoxDomain(const int64_t phaseSpace[DIMENSION],
                                 const double latticeConst, const double cutoffRadius);

    /**
     * In this method, each processor will be bound to a cartesian coordinate.
     *
     * It first divide the simulation box into N pieces(sub-box) (N is the count of all processors).
     * And each processor will be bound to a sub-box, and tagged with a cartesian coordinate(x,y,z).
     */
    void decomposition();

    void exchangeInter(atom *_atom);

    void borderInter(atom *_atom);

    void exchangeAtomfirst(atom *_atom);

    void exchangeAtom(atom *_atom);

    void sendrho(atom *_atom);

    void sendDfEmbed(atom *_atom);

    void sendForce(atom *_atom);

    /**
     * get lower bound of current sub-box at a dimension.
     * @param dimension
     * @return
     */
    double getSubBoxLowerBounding(int dimension) const;

    /**
     * get upper bound of current sub-box at some dimension specificed by {@var dimension}.
     * @param dimension
     * @return
     */
    double getUpperBoxLowerBounding(int dimension) const;

    /**
     * get ghost length at dimension
     * @param index
     * @return
     */
    double getGhostLength(int index) const;

private:

    /** local information for current simulation sub-box. **/
    /**
     * the count of processors at each dimension.
     */
    int _grid_size[DIMENSION] = {0};

    /**
     * The cartesian coordinate of the sub-box bound to this processor.
     */
    int _coords[DIMENSION];

    /**
     * the rank ids of contiguous processors in space.
     */
    int _rank_id_neighbours[DIMENSION][2];

    /**bound of local sub-box**/
    double _boundingBoxMin[DIMENSION]; // the lower bound of current sub-box.
    double _boundingBoxMax[DIMENSION]; // the upper bound of current sub-box.

    double _ghostLength[DIMENSION];  // ghost length, which equals to the cutoff radius.
    double _ghostBoundingBoxMin[DIMENSION]; // the lower ghost bound of current sub-box.
    double _ghostBoundingBoxMax[DIMENSION]; // the upper ghost bound of current sub-box.

    /** mpi data struct. **/
    MPI_Datatype _mpi_Particle_data;
    MPI_Datatype _mpi_latParticle_data;

    std::vector<std::vector<int> > sendlist;
    std::vector<std::vector<int> > recvlist;

    std::vector<std::vector<int> > intersendlist;
    std::vector<std::vector<int> > interrecvlist;

    /**
     * get the lower bound of current sub-box at dimension d.
     */
    double getBoundingBoxMin(int dimension);

    /**
     * get the upper bound of current sub-box at dimension d.
     */
    double getBoundingBoxMax(int dimension);
};

#endif // CRYSTAL_MD_DOMAIN_DECOMPOSITION_H
