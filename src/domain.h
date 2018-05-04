//
// Created by baihe back to 2016-12-22.
//

#ifndef CRYSTAL_MD_DOMAIN_DECOMPOSITION_H
#define CRYSTAL_MD_DOMAIN_DECOMPOSITION_H

#include <mpi.h>
#include <vector>

#include "atom.h"
#include "pre_config.h"

#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR

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
 *
 *
 * variable naming:
 * 1.variables for real length(measured length) and real bounding(measured bounding) of sub-box or global box
 * have a prefix "meas" or "_meas", with double type.
 *
 * 2.variables for cartesian coordinate of box decomposition have a prefix "grid_coord" or "_grid_coord", with int type;
 * and the grid count of box decomposition at each dimension have a prefix "grid_size" or "_grid_size", with int type;
 *
 * 3.variables for lattice coordinate have prefix "lattice_coord" or "_lattice_coord", with int type.
 *
 * 4.variables for lattice count in box or in global box have a prefix "lattice_size" or "_lattice_size", with int type.
 */

typedef int _type_lattice_size;

class DomainDecomposition {
public:

    /**
     * global information for the simulation box.
     */
    double _meas_global_length[DIMENSION]; // todo private.

    /**
     * the lower and upper bound coordinate for the global simulation box.
     */
    double _grid_coord_global_box_low[DIMENSION], _grid_coord_global_box_up[DIMENSION];

    DomainDecomposition(const int64_t *phaseSpace, const double latticeConst,
                        const double cutoffRadius);

    ~DomainDecomposition();

    /**
     * In this method, each processor will be bound to a cartesian coordinate.
     *
     * It first divide the simulation box into N pieces(sub-box) (N is the count of all processors).
     * And each processor will be bound to a sub-box, and tagged with a cartesian coordinate(x,y,z).
     */
    DomainDecomposition *decomposition();

    /**
     * set length of global simulation box.
     * and set upper and lower bound of global simulation box.
     * @param phaseSpace  phase space, the unit length of simulation box.
     * @param latticeConst lattice const
     */
    DomainDecomposition *createGlobalDomain();

    /**
     * set bound for current sub-box.
     * @param phaseSpace  phase space, the unit length of simulation box.
     * @param latticeConst lattice const
     */
    DomainDecomposition *createLocalBoxDomain();

    /**
     * get global measured length of the simulation box at dimension {@var d}.
     * @param d dimension
     * @return box length at dimension d.
     */
    inline double getMeasuredGlobalLength(int d) const {
        return _meas_global_length[d];
    }


    void exchangeInter(atom *_atom);

    void borderInter(atom *_atom);

    void exchangeAtomfirst(atom *_atom);

    void exchangeAtom(atom *_atom);

    void sendrho(atom *_atom);

    void sendDfEmbed(atom *_atom);

    void sendForce(atom *_atom);

    /**
     * get lower bound of current sub-box at a dimension.
     */
    double getMeasuredSubBoxLowerBounding(int dimension) const;

    /**
     * get upper bound of current sub-box at some dimension specified by {@var dimension}.
     */
    double getMeasuredSubBoxUpperBounding(int dimension) const;

    /**
     *  todo document
     */
    double getMeasuredGhostLowerBounding(int dimension) const;

    double getMeasuredGhostUpperBounding(int dimension) const;

    /**
     * get measured ghost length at dimension
     * @param index
     * @return
     */
    double getMeasuredGhostLength(int index) const;

    _type_lattice_size getSubBoxLatticeSize(unsigned short dimension) const;

private:

    const double _lattice_const, _cutoff_radius; // lattice const and cutoff radius for construct domain.
    int64_t _phase_space[DIMENSION]; // phase space of the global simulation box; the lattice count in each dimension in global simulation box.

    /** local information for current simulation sub-box. **/
    /**
     * the count of processors at each dimension.
     */
    int _grid_size[DIMENSION] = {0};

    /**
     * The cartesian coordinate of the sub-box bound to this processor after running grid decomposition.
     */
    int _grid_coord_sub_box[DIMENSION];

    /**
     * the rank ids of contiguous processors in space.
     */
    int _rank_id_neighbours[DIMENSION][2];

    /**bounding of local sub-box**/
    double _meas_sub_box_lower_bounding[DIMENSION]; // the measured lower bounding of current sub-box.
    double _meas_sub_box_upper_bounding[DIMENSION]; // the measured upper bounding of current sub-box.

    /**bounding of ghost of local sub-box**/
    double _meas_ghost_length[DIMENSION];  // measured ghost length, which equals to the cutoff radius.
    double _meas_ghost_lower_bounding[DIMENSION]; // the measured ghost lower bound of current sub-box.
    double _meas_ghost_upper_bounding[DIMENSION]; // the measured ghost upper bound of current sub-box.

    /*lattice count in local sub-box*/
    //  int _lattice_coord_lower[DIMENSION]; // the lower bounding of lattice coordinate for current sub-box.
    //  int _lattice_coord_upper[DIMENSION]; // the upper bounding of lattice coordinate for current sub-box.
    _type_lattice_size _lattice_size_sub_box[DIMENSION]; // lattice count of current sub-box in each dimension (upper bounding - lower bounding).

    /** mpi data struct. **/
    MPI_Datatype _mpi_Particle_data;
    MPI_Datatype _mpi_latParticle_data;

    std::vector<std::vector<int> > sendlist;
    std::vector<std::vector<int> > recvlist;

    std::vector<std::vector<int> > intersendlist;
    std::vector<std::vector<int> > interrecvlist;

};

#endif // CRYSTAL_MD_DOMAIN_DECOMPOSITION_H
