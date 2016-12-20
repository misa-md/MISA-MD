#include <stdlib.h>
#include <mpi.h>
#include <iostream>
#include <algorithm>

#include "InterpolationObject.h"

InterpolationObject::InterpolationObject(){
	values = NULL;
}

InterpolationObject::~InterpolationObject(){
	delete[] values;
	for(int i = 0; i < n + 1; i++)
		delete spline[i];
	delete[] spline;
}

void InterpolationObject::initInterpolationObject(int _n, double _x0, double dx, double* data)
{
   n = _n;
   values = new double[n+1];
   invDx = 1.0/dx;
   x0 = _x0;
   for (int ii = 0; ii < n; ++ii)
      values[ii+1] = data[ii];
}

void InterpolationObject::bcastInterpolationObject(int rank){
	struct
	{
		int n;
		double x0, invDx;
	} buf;

	if (rank == 0)
	{
		buf.n     = n;
		buf.x0    = x0;
		buf.invDx = invDx;
	}
    MPI_Bcast(&buf, sizeof(buf), MPI_BYTE, 0, MPI_COMM_WORLD);

	if (rank != 0)
	{
		n = buf.n;
		x0 = buf.x0;
		invDx = buf.invDx;
		values = new double[buf.n+1];
	}
	MPI_Bcast(values, n+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void InterpolationObject::interpolatefile(){
	spline = new double *[n+1];
	for(int i = 0; i < n+1; i++)
		spline[i] = new double[7];
	for (int m = 1; m <= n; m++) spline[m][6] = values[m];

	spline[1][5] = spline[2][6] - spline[1][6];
	spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
	spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
	spline[n][5] = spline[n][6] - spline[n-1][6];

	for (int m = 3; m <= n-2; m++)
		spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) +
		8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;

	for (int m = 1; m <= n-1; m++) {
		spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) -
			2.0*spline[m][5] - spline[m+1][5];
		spline[m][3] = spline[m][5] + spline[m+1][5] -
			2.0*(spline[m+1][6]-spline[m][6]);
	}

	spline[n][4] = 0.0;
	spline[n][3] = 0.0;

	for (int m = 1; m <= n; m++) {
		spline[m][2] = spline[m][5] * invDx;
		spline[m][1] = 2.0*spline[m][4] * invDx;
		spline[m][0] = 3.0*spline[m][3] * invDx;
	}
}
