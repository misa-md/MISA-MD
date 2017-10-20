#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_
#include "atom.h"

class integrator{
public:
	integrator(double timestepLength);
	~integrator();

	void firststep(atom* _atom);
	void secondstep(atom* _atom);
	void setTimestepLength(double dt) {
		_timestepLength = dt;
	}
private:
	double _timestepLength;
};

#endif /* INTEGRATOR_H_ */
