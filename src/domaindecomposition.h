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
	void setGridSize(int num_procs);

	MPI_Comm _comm;
	int _comm_size;

	MPI_Datatype _mpi_Particle_data;
	MPI_Datatype _mpi_latParticle_data;
	//! 每个维度的进程数
	int _gridSize[DIM];
	//! 进程的笛卡尔坐标
	int _coords[DIM];
	//!  进程号
	int _rank;
	//! 邻居进程号
	int _neighbours[DIM][2];

	vector<vector<int> > sendlist;
	vector<vector<int> > recvlist;

	vector<vector<int> > intersendlist;
	vector<vector<int> > interrecvlist;
};

#endif /* DOMAINDECOMPOSITION_H_ */
