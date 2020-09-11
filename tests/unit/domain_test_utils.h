//
// Created by genshen on 5/12/18.
//

#ifndef MISA_MD_DOMAIN_TEST_UTILS_H
#define MISA_MD_DOMAIN_TEST_UTILS_H

#include <cstdint>
#include <comm/domain/bcc_domain.h>

comm::BccDomain *getDomainInstance(int64_t space[3], double lattice_const, double cutoff_radius);

#endif //MISA_MD_DOMAIN_TEST_UTILS_H
