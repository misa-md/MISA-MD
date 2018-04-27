//
// Created by baihe back to 2015-06-30.
//

#ifndef CRYSTAL_MD_DOMAIN_H
#define CRYSTAL_MD_DOMAIN_H

#include <utils/data_def.h>

class domain {

public:
    domain(kiwi::RID rank);

    double getGlobalLength(int index) const;

    void setGlobalLength(int index, double length);

    double getGlobalLength(int d) { return _globalLength[d]; }

private:
    int _localRank;
    double _globalLength[3];

    domain();

    domain(domain &domain);

    domain &operator=(domain &domain);
};

#endif //CRYSTAL_MD_DOMAIN_H
