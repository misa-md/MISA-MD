#ifndef CREATEATOM_H
#define CREATEATOM_H

#include "atom.h"

#define boltz 8.617343e-5
#define mvv2e 1.0364269e-4

class createatom{
public:
	createatom(double tset);
	~createatom();
	void createphasespace(atom* _atom, double mass, int box_x, int box_y, int box_z);
private:
	double dof_compute(unsigned long natom);
	double t_set;
};

#endif
