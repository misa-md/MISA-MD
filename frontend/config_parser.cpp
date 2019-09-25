//
// Created by gensh(genshenchu@gmail.com) on 2017/4/16.
//

#include <iostream>
#include <fstream>
#include <utils/mpi_utils.h>
#include <utils/bundle.h>
#include <climits>
#include <cassert>

#include "config_parser.h"
#include "utils/rpcc.hpp"
#include "def_config_values.h"


ConfigParser *ConfigParser::m_pInstance = nullptr;

ConfigParser::ConfigParser() : kiwi::config::config() {
    configValues = ConfigValues(); // todo ConfigValues cap
}

// a simple single mode.
ConfigParser *ConfigParser::getInstance() {
    if (m_pInstance == nullptr) {
        m_pInstance = new ConfigParser();
    }
    return m_pInstance; // make sure there is a configure instance.
}

ConfigParser *ConfigParser::newInstance(const std::string &configureFilePath) {
    if (m_pInstance == nullptr) {
        m_pInstance = new ConfigParser();  // todo delete
        m_pInstance->parseConfig(configureFilePath);
    }
    return m_pInstance;
}

bool ConfigParser::configureCheck() {
    //todo
    return true;
}

// @override only for master processor.
// @see https://github.com/jbeder/yaml-cpp/wiki for more details.
void ConfigParser::parseConfig(const std::string config_file) {
    std::ifstream ifs(config_file);
    if (!ifs.good()) { // todo fixme important, if file not exist, this branch would not enter.
        setError("can not access the configure file: " + config_file);
        return;
    }
    const YAML::Node config = YAML::Load(ifs);
    ifs.close();

    // parse simulation
    if (!parseConfigSimulation(config["simulation"])) {
        return;
    }

    // parse creation
    if (!parseConfigCreation(config["creation"])) {
        return;
    }

    // parse potential
    if (!parseConfigPotential(config["potential"])) {
        return;
    }

    // parse output
    if (!parseConfigOutput(config["output"])) {
        return;
    }

    // parse stages
    if (!parseStages(config["stages"])) {
        return;
    }
}

// @override
void ConfigParser::putConfigData(kiwi::Bundle &bundle) {
    configValues.packdata(bundle);
}

// @override
void ConfigParser::getConfigData(kiwi::Bundle &bundle) {
    configValues.unpackdata(bundle);
}

bool ConfigParser::parseConfigSimulation(const YAML::Node &yaml_sim) {
    if (!yaml_sim) {
        setError("\"simulation\" in config is not specified.");
        return false;
    }
    // resolve simulation.phase_space
    if (!yaml_sim["phasespace"].IsSequence() || yaml_sim["phasespace"].size() != DIMENSION) {
        setError("\"simulation.phasespace\" must be a array with length 3.");
        return false;
    } else {
        for (int i = 0; i < DIMENSION; i++) {
            configValues.phaseSpace[i] = yaml_sim["phasespace"][i].as<unsigned long>(0); // todo uint64
        }
    }

    //resolve simulation.cutoff_radius
    if (!yaml_sim["cutoff_radius_factor"]) {
        setError("cutoff_radius_factor must be set in config.");
        return false;
    } else {
        configValues.cutoffRadiusFactor = yaml_sim["cutoff_radius_factor"].as<double>(0.0);
    }

    //resolve simulation.latticeconst
    if (!yaml_sim["lattice_const"]) {
        setError("lattice const must be set in config.");
        return false;
    } else {
        configValues.latticeConst = yaml_sim["lattice_const"].as<double>(0.0);
    }

    configValues.timeStepLength = yaml_sim["def_timesteps_length"].as<double>(default_time_length);
    return true;
}

bool ConfigParser::parseConfigCreation(const YAML::Node &yaml_creation) {
    if (!yaml_creation) {
        setError("\"creation\" in config is not specified.");
        return false;
    }

    configValues.createPhaseMode = yaml_creation["create_phase"].as<bool>(default_create_phase);// default value is true
    if (configValues.createPhaseMode) { //create mode
        if (yaml_creation["create_t_set"]) {
            configValues.createTSet = yaml_creation["create_t_set"].as<double>(0);
        } else {
            setError("creation t_set must be specified.");
            return false;
        }
        configValues.createSeed = yaml_creation["create_seed"].as<int>(default_random_seek);
    } else {  // read mode.
        auto yaml_read_file = yaml_creation["read_phase_filename"];
        if (yaml_read_file) {
            configValues.readPhaseFilename = yaml_read_file.as<std::string>();
        } else {
            setError("read phase file must be specified.");
            return false;
        }
    }

    // resolve simulation.alloy
    return parseConfigAlloy(yaml_creation["alloy"]);
}

bool ConfigParser::parseConfigPotential(const YAML::Node &yaml_pot) {
    //potential_file
    if (!yaml_pot) {
        setError("\"potential\" in config is not specified.");
        return false;
    }
    const YAML::Node yaml_pot_type = yaml_pot["type"];
    const YAML::Node yaml_pot_file = yaml_pot["file_path"];

    if (yaml_pot_type) {
        configValues.potentialFileType = yaml_pot_type.as<std::string>();
    } else {
        setError("potential file type must be specified.");
        return false;
    }

    if (yaml_pot_file) {
        configValues.potentialFilename = yaml_pot_file.as<std::string>();;
    } else {
        setError("potential filename must be specified.");
        return false;
    }
    return true;
}

bool ConfigParser::parseConfigOutput(const YAML::Node &yaml_output) {
    if (!yaml_output) {
        setError("output section in config file is not specified.");
        return false;
    }
    // resole dump config.
    const YAML::Node yaml_dump = yaml_output["dump"];
    if (!yaml_dump) {
        setError("dumping output section in config file is not specified.");
        return false;
    }

    // dump mode
    if (yaml_dump["atoms_dump_mode"]) {
        if ("copy" == yaml_dump["atoms_dump_mode"].as<std::string>("")) { // todo equal?
            configValues.output.atomsDumpMode = OutputMode::COPY;
        } else {
            configValues.output.atomsDumpMode = OutputMode::DEBUG;
        }
    }
    // todo 0 as default.
    configValues.output.atomsDumpInterval = yaml_dump["atoms_dump_interval"].as<uint64_t>(ULONG_MAX);
    configValues.output.outByFrame = yaml_dump["by_frame"].as<bool>(false);

    configValues.output.atomsDumpFilePath = yaml_dump["atoms_dump_file_path"].as<std::string>(
            DEFAULT_OUTPUT_DUMP_FILE_PATH);
    configValues.output.originDumpPath = yaml_dump["origin_dump_path"].as<std::string>("");
    // todo check if it is a real path.
    if (configValues.output.outByFrame && configValues.output.atomsDumpFilePath.find("{}") == std::string::npos) {
        setError("error format of dump file path");
    }

    // resole thermodynamics config
    const YAML::Node yaml_thermo = yaml_output["thermo"];
    if (!yaml_thermo) {
        setError("thermodynamics output section in config file is not specified.");
        return false;
    }
    configValues.output.thermo_interval = yaml_thermo["interval"].as<uint64_t>(0);

    // resolve logs.
    const YAML::Node yaml_logs = yaml_output["logs"];
    if (!yaml_logs) {
        setError("output logs in config file is not specified.");
        return false;
    }

    std::string logs_mode_str = yaml_logs["logs_mode"].as<std::string>(DEFAULT_LOGS_MODE_CONSOLE_STRING);
    if (LOGS_MODE_CONSOLE_STRING == logs_mode_str) {
        configValues.output.logs_mode = LOGS_MODE_CONSOLE;
    } else {
        configValues.output.logs_mode = LOGS_MODE_FILE;
        configValues.output.logs_filename = yaml_logs["logs_filename"].as<std::string>("");
        // if filename is empty, we will generate a filename.
        if (configValues.output.logs_filename.empty()) {
            std::ostringstream str_stream;
            str_stream << "md.";
            str_stream << rpcc();
            str_stream << ".log";
            configValues.output.logs_filename = str_stream.str();
        }
    }
    return true;
}

bool ConfigParser::parseConfigAlloy(const YAML::Node &yaml_alloy) {
    if (!yaml_alloy) {
        setError("alloy config is not specified in creation mode.");
        return false;
    }
    configValues.alloyCreateSeed = yaml_alloy["create_seed"].as<int>(default_random_seek);

    const YAML::Node ratios = yaml_alloy["ratio"];
    if (!ratios.IsMap()) {
        setError("alloy ratio must be a map.");
        return false;
    } else {
        configValues.alloyRatio[atom_type::Fe] = ratios["Fe"].as<int>(1);
        configValues.alloyRatio[atom_type::Cu] = ratios["Cu"].as<int>(1);
        configValues.alloyRatio[atom_type::Ni] = ratios["Nii"].as<int>(1);
    }
    return true;
}

bool ConfigParser::resolveConfigCollision(Stage *stage, const YAML::Node &yaml_collision) {
    if (!yaml_collision) {
        setError("collision config is not set.");
        return false;
    }
    stage->collisionStep = yaml_collision["collision_step"].as<unsigned long>(0); // todo default 0.

    const YAML::Node yaml_lat = yaml_collision["lat"];
    if (yaml_lat.IsSequence() && yaml_lat.size() == 4) {
        for (std::size_t i = 0; i < 4; i++) {
            stage->collisionLat[i] = yaml_lat[i].as<int64_t>();  // todo conversion.
        }
    } else { //the array length must be 4.
        setError("array length of value \"collision.lat\" must be 4.");
        return false;
    }
    stage->pkaEnergy = yaml_collision["energy"].as<double>(0.0);


    const YAML::Node yaml_dir = yaml_collision["direction"];
    if (yaml_dir.IsSequence() && yaml_dir.size() == DIMENSION) {
        for (std::size_t i = 0; i < DIMENSION; i++) {
            stage->direction[i] = yaml_dir[i].as<double>(0.0);  // todo conversion.
        }
    } else { //the array length must be 3.
        setError("array length of value \"simulation.collision.collision_v\" must be 3.");
        return false;
    }
    return true;
}

bool ConfigParser::parseStages(const YAML::Node &yaml_stages) {
    if (!yaml_stages && yaml_stages.IsSequence()) {
        setError("stages is not correctly set in config file.");
        return false;
    }

    //resolve simulation.collision
//    resolveConfigCollision(table->get_table("collision"));
    std::vector<unsigned long> vsl_break_points;
    std::vector<double> vsl_lengths;
    const std::size_t stages_size = yaml_stages.size();
    for (std::size_t i = 0; i < stages_size; i++) {
        const YAML::Node yaml_stage = yaml_stages[i];
        Stage stage;
        stage.steps = yaml_stage["steps"].as<unsigned long>(0);
        stage.step_length = yaml_stage["step_length"].as<double>(0.0);
        if (yaml_stage["set_v"]) {
            stage.collision_set = true;
            resolveConfigCollision(&stage, yaml_stage["set_v"]);
        } else {
            stage.collision_set = false;
        }
        configValues.stages.emplace_back(stage);
        configValues.timeSteps += stage.steps; // calculate total steps.
    }
    return true;
}
