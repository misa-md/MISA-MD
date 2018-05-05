#include "create_atom.h"

#include <mpi.h>
#include <cmath>

create_atom::create_atom(double tset) {
    t_set = tset;
}

create_atom::~create_atom() = default;

void create_atom::createphasespace(atom *_atom, double mass, int box_x, int box_y, int box_z) {
    double factor = 1.0 / sqrt(mass);
    _atom->createphasespace(factor, box_x, box_y, box_z);

    int nlocalatom = _atom->getnlocalatom();
    unsigned long natom = 2 * (unsigned long)box_x * (unsigned long)box_y * (unsigned long)box_z;
    double masstotal = natom * mass; // todo multiple type.
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

double create_atom::dof_compute(unsigned long natom) {
    double tfactor;
    unsigned long dof = 3 * natom;
    dof -= 3;
    tfactor = mvv2e / (dof * boltz);
    return tfactor;
}
