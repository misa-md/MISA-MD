#include "domain.h"

domain::domain(int rank){
	_localRank = rank;
	this->_globalLength[0] = 0;
	this->_globalLength[1] = 0;
	this->_globalLength[2] = 0;
}

void domain::setGlobalLength(int index, double length) {
	_globalLength[index] = length;
}
