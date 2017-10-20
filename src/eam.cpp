#include <mpi.h>
#include <cstring>
#include <stdio.h>

#include "eam.h"

using namespace std;

eam::eam(){};
eam::~eam(){
	delete[] f;
	delete[] phi;
	delete[] rho;
//	delete[] mass;
};

void eam::setatomicNo(double nAtomic){
	atomicNo = nAtomic;
}

void eam::setlat(double latticeconst){
	lat = latticeconst;
}

void eam::setmass(int i, double _mass){
	mass[i] = _mass;
}

void eam::setlatticeType(char* _latticeType){
	strcpy(latticeType, _latticeType);
}

void eam::setname(char* _name){
	strcpy(name, _name);
}

void eam::setcutoff(double _cutoff){
	cutoff = _cutoff;
}

void eam::init(int nElems){
	_nElems = nElems;
	f = new InterpolationObject[nElems];
	phi = new InterpolationObject[nElems + nElems * (nElems-1) / 2];
	rho = new InterpolationObject[nElems];
	mass = new double[nElems];
}

void eam::initf(int i, int nRho, double x0, double dRho, double* buf){
	f[i].initInterpolationObject(nRho, x0, dRho, buf);
}

void eam::initphi(int i, int nR, double x0, double dR, double *buf){
	phi[i].initInterpolationObject(nR, x0, dR, buf);
}

void eam::initrho(int i, int nR, double x0, double dR, double *buf){
	rho[i].initInterpolationObject(nR, x0, dR, buf);
}

void eam::eamBcast(int rank){
	MPI_Bcast(&_nElems, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank != 0){
		this->init(_nElems);
	}
	MPI_Bcast(&mass, _nElems, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	for (int i = 0; i < _nElems; i++){
		rho[i].bcastInterpolationObject(rank);
		f[i].bcastInterpolationObject(rank);
	}
	for (int i = 0; i< (_nElems + _nElems * (_nElems-1) / 2); i++)
		phi[i].bcastInterpolationObject(rank);
}

void eam::interpolatefile(){
	for (int i = 0; i < _nElems; i++){
		rho[i].interpolatefile();
		f[i].interpolatefile();
	}
	for (int i = 0; i < (_nElems + _nElems * (_nElems-1) / 2); i++)
		phi[i].interpolatefile();
}
