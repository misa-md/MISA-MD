#ifndef DOMAINDECOMPOSITION_H_
#define DOMAINDECOMPOSITION_H_

#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include <mpi.h>
#include <vector>

class atom;

#include "atom.h"
#include "domain.h"

using namespace std;

#define DIM 3

#define LOWER  0
#define HIGHER 1

class domaindecomposition{
public:
	/**
	 * In this construction method, each processor will be bound to a cartesian coordinate.
	 *
	 * It first divide the simulation space into N pieces(sub-space) (N is the count of all processors).
	 * First, if N can be decomposed as N = N_x * N_y * N_z, where N, N_x, N_y, N_z are all integer bigger than or equal to 1,
	 * then the whole simulation space will be divided into N sub-space with N_z levels in z axis,
	 * and in each level, it has  N_x * N_y sub-space.
	 *
	 * Second, we can bind each processor to a cartesian coordinate (x,y,z) due to the space partition, where 0 <= x < N_x, 0 <= z < N_z, 0 <= z < N_z.
	 */
	domaindecomposition();

	~domaindecomposition();

	void exchangeInter(atom* _atom, domain* domain);

	void borderInter(atom* _atom, domain* domain);

	void exchangeAtomfirst(atom* _atom, domain* domain);
	void exchangeAtom(atom* _atom, domain* domain);

	void sendrho(atom* _atom);
	void sendDfEmbed(atom* _atom);

	void sendforce(atom* _atom);

	double getBoundingBoxMin(int dimension, domain* domain);
	double getBoundingBoxMax(int dimension, domain* domain);
private:

	MPI_Datatype _mpi_Particle_data;
	MPI_Datatype _mpi_latParticle_data;
	// 每个维度的进程数
	int _gridSize[DIM];
	// 进程的笛卡尔坐标
	int _coords[DIM];
	// rank in communication _comm.
	int _rank;
	// 邻居进程号
	int _neighbours[DIM][2];

	vector<vector<int> > sendlist;
	vector<vector<int> > recvlist;

	vector<vector<int> > intersendlist;
	vector<vector<int> > interrecvlist;
};

#endif /* DOMAINDECOMPOSITION_H_ */
