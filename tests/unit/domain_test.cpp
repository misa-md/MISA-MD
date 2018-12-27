//
// Created by genshen on 2018-05-02.
//

#include <gtest/gtest.h>
#include <utils/mpi_utils.h>
#include <iostream>
#include <cmath>
#include "utils/mpi_domain.h"
#include "atom.h"
#include "domain_test_utils.h"

TEST(domain_test_decomposition, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    int rand_seek = 1024;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);
    auto *_atom = new atom(_domain, lattice_const, cutoff_radius_factor, rand_seek);

    int grid_sum[3] = {0, 0, 0};
    int local_grid[3] = {(_domain->lattice_size_sub_box[0]), (_domain->lattice_size_sub_box[1]),
                         (_domain->lattice_size_sub_box[2])};
    MPI_Reduce(local_grid, grid_sum, 3, MPI_INT, MPI_SUM, MASTER_PROCESSOR, MPIDomain::sim_processor.comm);

    if (MPIDomain::sim_processor.own_rank == MASTER_PROCESSOR) {
        EXPECT_EQ(grid_sum[1], (2 * space[1])); // todo test
    }

    delete _domain;
    delete _atom;
}

TEST(domain_local_lattice_size, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    int nlocalx = floor(_domain->meas_sub_box_upper_bounding[0] / (lattice_const)) -
                  floor(_domain->meas_sub_box_lower_bounding[0] / (lattice_const));
    nlocalx *= 2;
    int nlocaly = floor(_domain->meas_sub_box_upper_bounding[1] / lattice_const) -
                  floor(_domain->meas_sub_box_lower_bounding[1] / lattice_const);
    int nlocalz = floor(_domain->meas_sub_box_upper_bounding[2] / lattice_const) -
                  floor(_domain->meas_sub_box_lower_bounding[2] / lattice_const);

    EXPECT_EQ(nlocalx, {
        unsigned short dimension = 0;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_size_sub_box[dimension];
    });
    EXPECT_EQ(nlocaly, {
        unsigned short dimension = 1;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_size_sub_box[dimension];
    });
    EXPECT_EQ(nlocalz, {
        unsigned short dimension = 2;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_size_sub_box[dimension];
    });

    delete _domain;
}

TEST(domain_ghost_lattice_size, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    int nghostx = _domain->lattice_size_sub_box[0] + 2 * 2 * ceil(_domain->meas_ghost_length[0] / lattice_const);
    int nghosty = _domain->lattice_size_sub_box[1] + 2 * ceil(_domain->meas_ghost_length[1] / lattice_const);
    int nghostz = _domain->lattice_size_sub_box[2] + 2 * ceil(_domain->meas_ghost_length[2] / lattice_const);

    EXPECT_EQ(nghostx, {
        unsigned short dimension = 0;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_size_ghost_extended[dimension];
    });
    EXPECT_EQ(nghosty, {
        unsigned short dimension = 1;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_size_ghost_extended[dimension];
    });
    EXPECT_EQ(nghostz, {
        unsigned short dimension = 2;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_size_ghost_extended[dimension];
    });

    delete _domain;
}

TEST(domain_local_lattice_coord, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    // lower boundary of lattice coordinate of local sub-box
    int lolocalx = floor(_domain->meas_sub_box_lower_bounding[0] / lattice_const) * 2;
    int lolocaly = floor(_domain->meas_sub_box_lower_bounding[1] / lattice_const);
    int lolocalz = floor(_domain->meas_sub_box_lower_bounding[2] / lattice_const);

    EXPECT_EQ(lolocalx, {
        unsigned short dimension = 0;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_coord_sub_box_lower[dimension];
    });
    EXPECT_EQ(lolocaly, {
        unsigned short dimension = 1;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_coord_sub_box_lower[dimension];
    });
    EXPECT_EQ(lolocalz, {
        unsigned short dimension = 2;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_coord_sub_box_lower[dimension];
    });

    // upper boundary of lattice coordinate of local sub-box
    int uplocalx = floor(_domain->meas_sub_box_upper_bounding[0] / lattice_const) * 2;
    int uplocaly = floor(_domain->meas_sub_box_upper_bounding[1] / lattice_const);
    int uplocalz = floor(_domain->meas_sub_box_upper_bounding[2] / lattice_const);

    EXPECT_EQ(uplocalx, {
        unsigned short dimension = 0;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_coord_sub_box_upper[dimension];
    });
    EXPECT_EQ(uplocaly, {
        unsigned short dimension = 1;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_coord_sub_box_upper[dimension];
    });
    EXPECT_EQ(uplocalz, {
        unsigned short dimension = 2;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_coord_sub_box_upper[dimension];
    });

    delete _domain;
}

TEST(domain_ghost_lattice_coord, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    // lower boundary of lattice coordinate of ghost
    int loghostx = _domain->lattice_coord_sub_box_lower[0] - 2 * ceil(cutoff_radius_factor / lattice_const);
    int loghosty = _domain->lattice_coord_sub_box_lower[1] - ceil(cutoff_radius_factor / lattice_const);
    int loghostz = _domain->lattice_coord_sub_box_lower[2] - ceil(cutoff_radius_factor / lattice_const);

    EXPECT_EQ(loghostx, {
        unsigned short dimension = 0;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_coord_ghost_lower[dimension];
    });
    EXPECT_EQ(loghosty, {
        unsigned short dimension = 1;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_coord_ghost_lower[dimension];
    });
    EXPECT_EQ(loghostz, {
        unsigned short dimension = 2;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_coord_ghost_lower[dimension];
    });

    // upper boundary of lattice coordinate of ghost
    int upghostx = _domain->lattice_coord_sub_box_upper[0] + 2 * ceil(cutoff_radius_factor / lattice_const);
    int upghosty = _domain->lattice_coord_sub_box_upper[1] + ceil(cutoff_radius_factor / lattice_const);
    int upghostz = _domain->lattice_coord_sub_box_upper[2] + ceil(cutoff_radius_factor / lattice_const);

    EXPECT_EQ(upghostx, {
        unsigned short dimension = 0;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_coord_ghost_upper[dimension];
    });
    EXPECT_EQ(upghosty, {
        unsigned short dimension = 1;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_coord_ghost_upper[dimension];
    });
    EXPECT_EQ(upghostz, {
        unsigned short dimension = 2;
        Domain *receiver = _domain;
        _type_lattice_size result;
        result= receiver->lattice_coord_ghost_upper[dimension];
    });

    delete _domain;
}
