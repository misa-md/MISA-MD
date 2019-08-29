//
// Created by genshen on 2019-08-28.
//

#include <gtest/gtest.h>
#include <utils/mpi_utils.h>
#include "toml_config.h"
#include "../test_config.h"

const std::string vsl_test_config_header = R"(
[simulation]
phasespace = [50, 50, 50]
cutoff_radius_factor = 1.96125
lattice_const = 2.85532
)";

const std::string vsl_test_config_tail = R"(
[simulation.createphase]
    create_phase = true
    create_t_set = 600.0
    create_seed = 466953

    [simulation.alloy]
        create_seed = 1024
        [simulation.alloy.ratio]
            Fe = 97
            Cu = 2
            Ni = 1

    [simulation.collision]
    collision_step = 2
    lat = [2, 2, 2, 0]
    pka = 6.8
    direction = [1.0, 1.0, 1.0]

    [simulation.potential_file]
    type = "setfl"
    filename = "FeCuNi.eam.alloy"

[output]
atoms_dump_mode = "copy"
atoms_dump_file_path = "crystal_mdl.{}.out"
origin_dump_path = ""
atoms_dump_interval = 10
by_frame=true
    [output.logs]
    logs_mode = "console"
    logs_filename = ""
)";

// test variable step length in config parsing level.
// @MPI
TEST(config_parsing_vsl_test, config_parsing_test) {
    const std::string config_filePath = TEST_TEMP_FILES_DIR "/test_config_parsing_vsl.toml";
    std::fstream fs(config_filePath, std::ios::out);
    if (!fs.good()) {
        GTEST_FATAL_FAILURE_("open file fail.");
    }
    fs << vsl_test_config_header << R"(
timesteps =  10
def_timesteps_length = 0.001
variable_step_length = [ [1, 4, 8], [0.0001, 0.0005, 0.001] ]
)"
       << vsl_test_config_tail;
    fs.close();

    ConfigParser *p_config;
    if (kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) {
        // initial config Obj, then read and resolve config file.
        p_config = ConfigParser::newInstance(config_filePath); // todo config file from argv.
        if (p_config->hasError) {
            GTEST_FATAL_FAILURE_(p_config->errorMessage.c_str());
        }
    } else {
        // just initial a empty config Obj.
        p_config = ConfigParser::getInstance();
    }
    p_config->sync(); // sync config data to other processors from master processor.

    // test
    const int expected_vsl_size = 3;
    const unsigned long expected_breaks[] = {1, 4, 8};
    const double expected_vsl_length[] = {0.0001, 0.0005, 0.001};

    EXPECT_EQ(p_config->configValues.vsl_size, expected_vsl_size);
    for (int i = 0; i < expected_vsl_size; i++) {
        EXPECT_EQ(p_config->configValues.vsl_break_points[i], expected_breaks[i]);
        EXPECT_EQ(p_config->configValues.vsl_lengths[i], expected_vsl_length[i]);
    }
}

// test variable step length in config parsing level.
// if variable step length is not set.
// @MPI
TEST(config_parsing_non_vsl_test, config_parsing_test) {
    const std::string config_filePath = TEST_TEMP_FILES_DIR "/test_config_parsing_vsl.toml";
    std::fstream fs(config_filePath, std::ios::out);
    if (!fs.good()) {
        GTEST_FATAL_FAILURE_("open file fail.");
    }
    fs << vsl_test_config_header << R"(
timesteps =  10
def_timesteps_length = 0.001
)"
       << vsl_test_config_tail;
    fs.close();

    ConfigParser *p_config;
    if (kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) {
        // initial config Obj, then read and resolve config file.
        p_config = ConfigParser::newInstance(config_filePath); // todo config file from argv.
        if (p_config->hasError) {
            GTEST_FATAL_FAILURE_(p_config->errorMessage.c_str());
        }
    } else {
        // just initial a empty config Obj.
        p_config = ConfigParser::getInstance();
    }
    p_config->sync(); // sync config data to other processors from master processor.

    // do assert
    EXPECT_EQ(p_config->configValues.vsl_size, 0);
    EXPECT_EQ(p_config->configValues.vsl_break_points.size(), 0);
    EXPECT_EQ(p_config->configValues.vsl_lengths.size(), 0);
}