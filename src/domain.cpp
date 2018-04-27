#include "domain.h"

domain::domain(kiwi::RID rank) {
    _localRank = rank;
    _globalLength[0] = 0;
    _globalLength[1] = 0;
    _globalLength[2] = 0;
}

void domain::setGlobalLength(int index, double length) {
    _globalLength[index] = length;
}
