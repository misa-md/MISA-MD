//
// Created by genshen on 5/12/18.
//

#ifndef MISA_MD_DOMAIN_TEST_FIXTURE_H
#define MISA_MD_DOMAIN_TEST_FIXTURE_H

#include <cstdint>
#include <comm/domain/bcc_domain.h>
#include <gtest/gtest.h>

class DomainFixture : public ::testing::Test {
protected:
    void SetUp() override;

private:
    static comm::BccDomain *getDomainInstance(int64_t space[3], double lattice_const, double cutoff_radius);

protected:
    const int64_t space[3] = {50, 60, 72};
    const double lattice_const = 0.86;
    const double cutoff_radius_factor = 1.1421;

    comm::BccDomain *p_domain = nullptr;
};

#endif //MISA_MD_DOMAIN_TEST_FIXTURE_H
