#include "createatom.h"

#include <mpi.h>
#include <cmath>

createatom::createatom(double tset) {
    t_set = tset;
}

createatom::~createatom() {}

void createatom::createphasespace(atom *_atom, double mass, int box_x, int box_y, int box_z) {
    double factor = 1.0 / sqrt(mass);
    _atom->createphasespace(factor, box_x, box_y, box_z);

    int nlocalatom = _atom->getnlocalatom();
    int natom = box_x * 2 * box_y * box_z;
    double masstotal = natom * mass;
    double p[3], vcm[3];
    p[0] = p[1] = p[2] = 0.0;
    _atom->vcm(mass, masstotal, p);
    MPI_Allreduce(p, vcm, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (masstotal > 0.0) {
        vcm[0] /= masstotal;
        vcm[1] /= masstotal;
        vcm[2] /= masstotal;
    }

    _atom->zero_momentum(vcm);

    double t = _atom->compute_scalar(mass);
    double scalar, tfactor;
    //double t_test, scalar_test;
    //t_test = t;
    MPI_Allreduce(&t, &scalar, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    tfactor = dof_compute(natom);
    scalar *= tfactor;
    //t_test *= tfactor;
    //MPI_Allreduce(&t_test,&scalar_test,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    _atom->rescale(scalar, t_set);
}

double createatom::dof_compute(unsigned long natom) {
    double tfactor;
    unsigned long dof = 3 * natom;
    dof -= 3;
    tfactor = mvv2e / (dof * boltz);
    return tfactor;
}
