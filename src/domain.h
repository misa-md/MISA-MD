//
// Created by baihe back to 2016-12-22.
//

#ifndef CRYSTAL_MD_DOMAIN_DECOMPOSITION_H
#define CRYSTAL_MD_DOMAIN_DECOMPOSITION_H

#include <mpi.h>
#include <vector>

#include "atom.h"
#include "pre_define.h"

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
 * 1.variables for real length(measured length) and real boundary(measured boundary) of sub-box or global box
 * have a prefix "meas" or "_meas", with double type.
 *
 * 2.variables for cartesian coordinate of box decomposition have a prefix "grid_coord" or "_grid_coord", with int type;
 * and the grid count of box decomposition at each dimension have a prefix "grid_size" or "_grid_size", with int type;
 *
 * 3.variables for lattice coordinate have prefix "lattice_coord" or "_lattice_coord", with int type.
 *
 * 4.variables for lattice count in box or in global box have a prefix "lattice_size" or "_lattice_size", with int type.
 */

class Domain {
public:

    Domain(const int64_t *phaseSpace, const double latticeConst, const double cutoffRadius);

    ~Domain();

    /**
     * In this method, each processor will be bound to a cartesian coordinate.
     *
     * It first divide the simulation box into N pieces(sub-box) (N is the count of all processors).
     * And each processor will be bound to a sub-box, and tagged with a cartesian coordinate(x,y,z).
     */
    Domain *decomposition();

    /**
     * set length of global simulation box.
     * and set upper and lower bound of global simulation box.
     * @param phaseSpace  phase space, the unit length of simulation box.
     * @param latticeConst lattice const
     */
    Domain *createGlobalDomain();

    /**
     * set boundary for current sub-box.
     *
     * Note: why some variable in x dimension is multiplied by 2:
     * the basic unit in x dimension is half the lattice, not a lattice like in y,z dimension.
     * in x dimension:
     *    | :  |  :  |  :  |  :  |  :  |  ...
     * x: 0 1  2  3  4  5  6  7  8  9  10 ...
     * but, in y or z dimension:
     *    | :  |  :  |  :  |  :  |  :  | ...
     * y: 0    1     2     3     4     5  ...
     * in above figure, |  :  | represents a lattice length.
     *
     * @param phaseSpace  phase space, the unit length of simulation box.
     * @param latticeConst lattice const
     */
    Domain *createSubBoxDomain();

    void exchangeInter(atom *_atom);

    void borderInter(atom *_atom);

    void exchangeAtomFirst(atom *_atom);

    void exchangeAtom(atom *_atom);

    void sendrho(atom *_atom);

    void sendDfEmbed(atom *_atom);

    void sendForce(atom *_atom);

    /**
     * get global measured length of the simulation box at dimension {@var d}.
     * @param d dimension
     * @return global simulation box length at dimension d.
     */
    inline double getMeasuredGlobalLength(int d) const {
        return _meas_global_length[d];
    }

    /**
     * get decomposed grid size at each dimension {@var d}.
     */
    inline double getDecomposedGridSize(int d) const {
        return _grid_size[d];
    }

    /**
     * get coordinate of lower boundary of global simulation box.
     */
    inline double getMeasuredGlobalBoxCoordLower(int d) const {
        return _meas_global_box_coord_lower[d];
    }

    /**
     * get coordinate of upper boundary of global simulation box.
     */
    inline double getMeasuredGlobalBoxCoordUpper(int d) const {
        return _meas_global_box_coord_upper[d];
    }

    /**
     * get measured lower bound of current sub-box at a dimension.
     */
    inline double getMeasuredSubBoxLowerBounding(int dimension) const {
        return _meas_sub_box_lower_bounding[dimension];
    }

    /**
     * get measured upper boundary of current sub-box at some dimension specified by {@var dimension}.
     */
    inline double getMeasuredSubBoxUpperBounding(int dimension) const {
        return _meas_sub_box_upper_bounding[dimension];
    }

    /**
     *  get measured ghost lower boundary of current sub-box.
     */
    inline double getMeasuredGhostLowerBounding(int dimension) const {
        return _meas_ghost_lower_bounding[dimension];
    }

    /**
     *  get measured ghost upper boundary of current sub-box.
     */
    inline double getMeasuredGhostUpperBounding(int dimension) const {
        return _meas_ghost_upper_bounding[dimension];
    }

    /**
     * get measured ghost length at dimension
     */
    inline double getMeasuredGhostLength(int index) const {
        return _meas_ghost_length[index];
    }

    /**
     * get lattice count in current sub-box.
     */
    inline _type_lattice_size getSubBoxLatticeSize(unsigned short dimension) const {
        return _lattice_size_sub_box[dimension];
    }

    /**
     * get lattice count in ghost area plus local area of current sub-box.
     */
    inline _type_lattice_size getGhostExtLatticeSize(unsigned short dimension) const {
        return _lattice_size_ghost_extended[dimension];
    }

    /**
     * get lattice count in ghost area plus local area of current sub-box.
     * which  @var _lattice_size_ghost_extended[d] = @var _lattice_size_sub_box[d] + 2 * @var_ lattice_size_ghost[d];
     * and also equals to _lattice_coord_sub_box_lower[d] - _lattice_coord_ghost_lower[d]
     */
    inline _type_lattice_size getGhostLatticeSize(unsigned short dimension) const {
        return _lattice_size_ghost[dimension];
    }

    /**
     * get lower boundary of lattice coordinate of current sub-box, in global coordinate system(GCY).
     */
    inline _type_lattice_size getGlobalSubBoxLatticeCoordLower(unsigned short dimension) const {
        return _lattice_coord_sub_box_lower[dimension];
    }

    /**
     * get upper boundary of lattice coordinate of current sub-box, in global coordinate system(GCY).
     */
    inline _type_lattice_size getGlobalSubBoxLatticeCoordUpper(unsigned short dimension) const {
        return _lattice_coord_sub_box_upper[dimension];
    }

    /**
     * get lower boundary of lattice coordinate of ghost of current sub-box, in global coordinate system(GCY).
     */
    inline _type_lattice_size getGlobalGhostLatticeCoordLower(unsigned short dimension) const {
        return _lattice_coord_ghost_lower[dimension];
    }

    /**
     * get upper boundary of lattice coordinate of ghost of current sub-box, in global coordinate system(GCY).
     */
    inline _type_lattice_size getGlobalGhostLatticeCoordUpper(unsigned short dimension) const {
        return _lattice_coord_ghost_upper[dimension];
    }

    /**
     * get lower boundary of lattice coordinate of current sub-box, in local coordinate system(LCY).
     */
    inline _type_lattice_size getLocalSubBoxLatticeCoordLower(unsigned short dimension) const {
        return _local_lattice_coord_sub_box_lower[dimension];
    }

    /**
     * get upper boundary of lattice coordinate of current sub-box, in local coordinate system(LCY).
     */
    inline _type_lattice_size getLocalSubBoxLatticeCoordUpper(unsigned short dimension) const {
        return _local_lattice_coord_sub_box_upper[dimension];
    }

    /**
     * get lower boundary of lattice coordinate of ghost of current sub-box, in local coordinate system(LCY).
     */
    inline _type_lattice_size getLocalGhostLatticeCoordLower(unsigned short dimension) const {
        return _local_lattice_coord_ghost_lower[dimension];
    }

    /**
     * get upper boundary of lattice coordinate of ghost of current sub-box, in local coordinate system(LCY).
     */
    inline _type_lattice_size getLocalGhostLatticeCoordUpper(unsigned short dimension) const {
        return _local_lattice_coord_ghost_upper[dimension];
    }

private:

    const double _lattice_const, _cutoff_radius_factor; // lattice const and cutoff radius for construct domain.
    int64_t _phase_space[DIMENSION]; // phase space of the global simulation box; the lattice count in each dimension in global simulation box.

    /**
    * global information for the simulation box.
    */
    double _meas_global_length[DIMENSION];

    /**
     * the count of processors at each dimension.
     */
    int _grid_size[DIMENSION] = {0};

    /**
     * the lower and upper boundary coordinate for the global simulation box.
     */
    double _meas_global_box_coord_lower[DIMENSION], _meas_global_box_coord_upper[DIMENSION];

    /** local information for current simulation sub-box. **/
    /**
     * The cartesian coordinate of the sub-box bound to this processor after running grid decomposition.
     */
    int _grid_coord_sub_box[DIMENSION];

    /**
     * the rank ids of contiguous processors in space.
     */
    int _rank_id_neighbours[DIMENSION][2];

    /**boundary of local sub-box**/
    // the measured lower and upper boundary of current sub-box (global box Coordinate System).
    double _meas_sub_box_lower_bounding[DIMENSION], _meas_sub_box_upper_bounding[DIMENSION];

    /**boundary of ghost of local sub-box**/
    double _meas_ghost_length[DIMENSION];  // measured ghost length, which equals to the cutoff radius.
    // the measured ghost lower and upper bound of current sub-box.
    double _meas_ghost_lower_bounding[DIMENSION], _meas_ghost_upper_bounding[DIMENSION];

    /*lattice count in local sub-box*/
    //  int _lattice_coord_lower[DIMENSION]; // the lower boundary of lattice coordinate for current sub-box.
    //  int _lattice_coord_upper[DIMENSION]; // the upper boundary of lattice coordinate for current sub-box.
    // lattice count in local sub-box area at each dimension (upper boundary - lower boundary).
    _type_lattice_size _lattice_size_sub_box[DIMENSION];
    // lattice count in ghost area plus sub-box area at each dimension (upper boundary - lower boundary).
    _type_lattice_size _lattice_size_ghost_extended[DIMENSION];
    // purge ghost size, just lattice count in ghost area.
    _type_lattice_size _lattice_size_ghost[DIMENSION];

    /*lattice boundary of local sub-box and ghost, but, the Coordinate System is still the global box.*/
    // lower and upper boundary of lattice coordinate of local sub-box at each dimension.
    _type_lattice_coord _lattice_coord_sub_box_lower[DIMENSION], _lattice_coord_sub_box_upper[DIMENSION];

    // lower and upper boundary of lattice coordinate in ghost area of current sub-box area at each dimension.
    _type_lattice_coord _lattice_coord_ghost_lower[DIMENSION], _lattice_coord_ghost_upper[DIMENSION];


    /*
     * lattice boundary of local sub-box and ghost, this is in local box Coordinate System(not global box.).
     * For convenience usage for atom index in local sub-box.
     *  | 0:ghost_lower                 | sub_box_lower                   | sub_box_upper               | ghost_upper
     *  |-------------------------------|---------...---------------------|-----------------------------|
     */
    // lower and upper boundary of lattice coordinate of local sub-box at each dimension in local coordinate system.
    _type_lattice_coord _local_lattice_coord_sub_box_lower[DIMENSION], _local_lattice_coord_sub_box_upper[DIMENSION];

    // lower and upper boundary of lattice coordinate in ghost area of current sub-box area at each dimension in local coordinate system.
    _type_lattice_coord _local_lattice_coord_ghost_lower[DIMENSION], _local_lattice_coord_ghost_upper[DIMENSION];

    /** mpi data struct. **/
    MPI_Datatype _mpi_Particle_data;
    MPI_Datatype _mpi_latParticle_data;

    std::vector<std::vector<int> > sendlist;
    std::vector<std::vector<int> > recvlist;

    std::vector<std::vector<int> > intersendlist;
    std::vector<std::vector<int> > interrecvlist;

    /**
     * set lattice coordinate boundary of current sub-box in global coordinate system(GCY).
     */
    void setSubBoxDomainGCS();

    /**
     * set lattice coordinate boundary of current sub-box in local coordinate system(LCY).
     */
    void setSubBoxDomainLCS(); // todo test.
};

#endif // CRYSTAL_MD_DOMAIN_DECOMPOSITION_H
