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

    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);
    auto *_atom = new atom(_domain, lattice_const, cutoff_radius_factor);

    int grid_sum[3] = {0, 0, 0};
    int local_grid[3] = {(_domain->lattice_size_sub_box[0]), (_domain->lattice_size_sub_box[1]),
                         (_domain->lattice_size_sub_box[2])};
    MPI_Reduce(local_grid, grid_sum, 3, MPI_INT, MPI_SUM, MASTER_PROCESSOR, kiwi::mpiUtils::global_process.comm);

    grid_sum[0] /= kiwi::mpiUtils::global_process.all_ranks / _domain->grid_size[0];
    grid_sum[1] /= kiwi::mpiUtils::global_process.all_ranks / _domain->grid_size[1];
    grid_sum[2] /= kiwi::mpiUtils::global_process.all_ranks / _domain->grid_size[2];
    if (kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) {
        EXPECT_EQ(grid_sum[0], space[0]);
        EXPECT_EQ(grid_sum[1], space[1]);
        EXPECT_EQ(grid_sum[2], space[2]);
    }

    delete _domain;
    delete _atom;
}

TEST(domain_local_lattice_size, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    int nlocalx = floor(_domain->meas_sub_box_region.x_high / (lattice_const)) -
                  floor(_domain->meas_sub_box_region.x_low / (lattice_const));
    int nlocaly = floor(_domain->meas_sub_box_region.y_high / lattice_const) -
                  floor(_domain->meas_sub_box_region.y_low / lattice_const);
    int nlocalz = floor(_domain->meas_sub_box_region.z_high / lattice_const) -
                  floor(_domain->meas_sub_box_region.z_low / lattice_const);

    EXPECT_EQ(nlocalx, _domain->lattice_size_sub_box[0]);
    EXPECT_EQ(nlocaly, _domain->lattice_size_sub_box[1]);
    EXPECT_EQ(nlocalz, _domain->lattice_size_sub_box[2]);

    delete _domain;
}

TEST(domain_ghost_lattice_size, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);
    // ceil(_domain->meas_ghost_length[0] / lattice_const) equals to ceil(_cutoff_radius_factor)
    EXPECT_EQ(_domain->lattice_size_ghost[0], ceil(cutoff_radius_factor));
    EXPECT_EQ(_domain->lattice_size_ghost[1], ceil(cutoff_radius_factor));
    EXPECT_EQ(_domain->lattice_size_ghost[2], ceil(cutoff_radius_factor));

    int nghostx = _domain->lattice_size_sub_box[0] + 2 * ceil(_domain->meas_ghost_length[0] / lattice_const);
    int nghosty = _domain->lattice_size_sub_box[1] + 2 * ceil(_domain->meas_ghost_length[1] / lattice_const);
    int nghostz = _domain->lattice_size_sub_box[2] + 2 * ceil(_domain->meas_ghost_length[2] / lattice_const);

    EXPECT_EQ(nghostx, _domain->lattice_size_ghost_extended[0]);
    EXPECT_EQ(nghosty, _domain->lattice_size_ghost_extended[1]);
    EXPECT_EQ(nghostz, _domain->lattice_size_ghost_extended[2]);

    delete _domain;
}

TEST(domain_local_lattice_coord, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);
    // test delta values
    EXPECT_EQ(_domain->lattice_size_sub_box[0],
              _domain->local_sub_box_lattice_coord_region.x_high - _domain->local_sub_box_lattice_coord_region.x_low);
    EXPECT_EQ(_domain->lattice_size_sub_box[1],
              _domain->local_sub_box_lattice_coord_region.y_high - _domain->local_sub_box_lattice_coord_region.y_low);
    EXPECT_EQ(_domain->lattice_size_sub_box[2],
              _domain->local_sub_box_lattice_coord_region.z_high - _domain->local_sub_box_lattice_coord_region.z_low);

    EXPECT_EQ(_domain->lattice_size_ghost_extended[0],
              _domain->local_ghost_lattice_coord_region.x_high - _domain->local_ghost_lattice_coord_region.x_low);
    EXPECT_EQ(_domain->lattice_size_ghost_extended[1],
              _domain->local_ghost_lattice_coord_region.y_high - _domain->local_ghost_lattice_coord_region.y_low);
    EXPECT_EQ(_domain->lattice_size_ghost_extended[2],
              _domain->local_ghost_lattice_coord_region.z_high - _domain->local_ghost_lattice_coord_region.z_low);

}

TEST(domain_global_lattice_coord, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    // lower boundary of lattice coordinate of local sub-box
    int lolocalx = floor(_domain->meas_sub_box_region.x_low / lattice_const);
    int lolocaly = floor(_domain->meas_sub_box_region.y_low / lattice_const);
    int lolocalz = floor(_domain->meas_sub_box_region.z_low / lattice_const);

    EXPECT_EQ(lolocalx, _domain->lattice_coord_sub_box_region.x_low);
    EXPECT_EQ(lolocaly, _domain->lattice_coord_sub_box_region.y_low);
    EXPECT_EQ(lolocalz, _domain->lattice_coord_sub_box_region.z_low);

    // upper boundary of lattice coordinate of local sub-box
    int uplocalx = floor(_domain->meas_sub_box_region.x_high / lattice_const);
    int uplocaly = floor(_domain->meas_sub_box_region.y_high / lattice_const);
    int uplocalz = floor(_domain->meas_sub_box_region.z_high / lattice_const);

    EXPECT_EQ(uplocalx, _domain->lattice_coord_sub_box_region.x_high);
    EXPECT_EQ(uplocaly, _domain->lattice_coord_sub_box_region.y_high);
    EXPECT_EQ(uplocalz, _domain->lattice_coord_sub_box_region.z_high);

    // test delta values
    EXPECT_EQ(_domain->lattice_size_sub_box[0],
              _domain->lattice_coord_sub_box_region.x_high - _domain->lattice_coord_sub_box_region.x_low);
    EXPECT_EQ(_domain->lattice_size_sub_box[1],
              _domain->lattice_coord_sub_box_region.y_high - _domain->lattice_coord_sub_box_region.y_low);
    EXPECT_EQ(_domain->lattice_size_sub_box[2],
              _domain->lattice_coord_sub_box_region.z_high - _domain->lattice_coord_sub_box_region.z_low);
    delete _domain;
}

TEST(domain_ghost_lattice_coord, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    // lower boundary of lattice coordinate of ghost
    int loghostx = _domain->lattice_coord_sub_box_region.x_low - ceil(cutoff_radius_factor / lattice_const);
    int loghosty = _domain->lattice_coord_sub_box_region.y_low - ceil(cutoff_radius_factor / lattice_const);
    int loghostz = _domain->lattice_coord_sub_box_region.z_low - ceil(cutoff_radius_factor / lattice_const);

    EXPECT_EQ(loghostx, _domain->lattice_coord_ghost_region.x_low);
    EXPECT_EQ(loghosty, _domain->lattice_coord_ghost_region.y_low);
    EXPECT_EQ(loghostz, _domain->lattice_coord_ghost_region.z_low);

    // upper boundary of lattice coordinate of ghost
    int upghostx = _domain->lattice_coord_sub_box_region.x_high + ceil(cutoff_radius_factor / lattice_const);
    int upghosty = _domain->lattice_coord_sub_box_region.y_high + ceil(cutoff_radius_factor / lattice_const);
    int upghostz = _domain->lattice_coord_sub_box_region.z_high + ceil(cutoff_radius_factor / lattice_const);

    EXPECT_EQ(upghostx, _domain->lattice_coord_ghost_region.x_high);
    EXPECT_EQ(upghosty, _domain->lattice_coord_ghost_region.y_high);
    EXPECT_EQ(upghostz, _domain->lattice_coord_ghost_region.z_high);

    // test delta values
    EXPECT_EQ(_domain->lattice_size_ghost_extended[0],
              _domain->lattice_coord_ghost_region.x_high - _domain->lattice_coord_ghost_region.x_low);
    EXPECT_EQ(_domain->lattice_size_ghost_extended[1],
              _domain->lattice_coord_ghost_region.y_high - _domain->lattice_coord_ghost_region.y_low);
    EXPECT_EQ(_domain->lattice_size_ghost_extended[2],
              _domain->lattice_coord_ghost_region.z_high - _domain->lattice_coord_ghost_region.z_low);
    delete _domain;
}
