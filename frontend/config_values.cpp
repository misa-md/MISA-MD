//
// Created by genshen(genshenchu@gmail.com) on 2018-3-9.
//

#include <iostream>
#include <cassert>
#include <mpi.h>
#include <utils/bundle.h>
#include "config_values.h"
#include "def_config_values.h"

Stage::Stage() : steps(0), step_length(default_time_length),
                 collision_set(false), collisionStep(0), collisionLat{0, 0, 0, 0},
                 pkaEnergy(0.0), direction{1.0, 1.0, 1.0},
                 velocity_set(false), velocity_step(0),
                 velocity_region{0, 0, 0, 0, 0, 0}, velocity_value{0.0, 0.0, 0.0},
                 rescales_set(false), rescale_t(0.0), rescale_every(1) {}

void Stage::packdata(kiwi::Bundle &bundle) {
    bundle.put(steps);
    bundle.put(step_length);
    // output dump
    bundle.put(dump_set);
    bundle.put(dump_preset_use);
    bundle.put(dump_every_steps);
    // collision
    bundle.put(collision_set);
    bundle.put(collisionStep);
    bundle.put(4, collisionLat);
    bundle.put(pkaEnergy);
    bundle.put(DIMENSION, direction);
    // velocity
    bundle.put(velocity_set);
    bundle.put(velocity_step);
    bundle.put(3, velocity_value);
    bundle.put(6, velocity_region);
    // rescale
    bundle.put(rescales_set);
    bundle.put(rescale_t);
    bundle.put(rescale_every);
}

void Stage::unnpackdata(int &cursor, kiwi::Bundle &bundle) {
    bundle.get(cursor, steps); // get size,
    bundle.get(cursor, step_length);
    // output dump
    bundle.get(cursor, dump_set);
    bundle.get(cursor, dump_preset_use);
    bundle.get(cursor, dump_every_steps);
    // collision
    bundle.get(cursor, collision_set);
    bundle.get(cursor, collisionStep);
    bundle.get(cursor, 4, collisionLat);
    bundle.get(cursor, pkaEnergy);
    bundle.get(cursor, DIMENSION, direction);
    // velocity
    bundle.get(cursor, velocity_set);
    bundle.get(cursor, velocity_step);
    bundle.get(cursor, 3, velocity_value);
    bundle.get(cursor, 6, velocity_region);
    // rescale
    bundle.get(cursor, rescales_set);
    bundle.get(cursor, rescale_t);
    bundle.get(cursor, rescale_every);
}

void DumpConfig::packdata(kiwi::Bundle &bundle) {
    bundle.put(name);
    bundle.put(2 * DIMENSION, region);
    bundle.put(mode);
    bundle.put(steps);
    bundle.put(file_path);
    bundle.put(by_frame);
    bundle.put(dump_mask);
    bundle.put(dump_whole_system);
}

void DumpConfig::unnpackdata(int &cursor, kiwi::Bundle &bundle) {
    bundle.get(cursor, name);
    bundle.get(cursor, 2 * DIMENSION, region);
    bundle.get(cursor, mode);
    bundle.get(cursor, steps);
    bundle.get(cursor, file_path);
    bundle.get(cursor, by_frame);
    bundle.get(cursor, dump_mask);
    bundle.get(cursor, dump_whole_system);
}


void AtomType::packdata(kiwi::Bundle &bundle) const {
    bundle.put(name);
    bundle.put(mass);
    bundle.put(weight);
}

void AtomType::unnpackdata(int &cursor, kiwi::Bundle &bundle) {
    bundle.get(cursor, name);
    bundle.get(cursor, mass);
    bundle.get(cursor, weight);
}

ConfigValues::ConfigValues() :
        phaseSpace{0, 0, 0}, cutoffRadiusFactor(0.0), latticeConst(0.0),
        timeSteps(0),
        createPhaseMode(true), createTSet(0.0), createSeed(1024), readPhaseFilename(""),
        alloyCreateSeed(1024), types(),
        output(), stages() {}
// todo potential type and filename initialize.

void ConfigValues::packdata(kiwi::Bundle &bundle) {
    // append data into buffer.
    bundle.put(DIMENSION, phaseSpace); // todo remove MPI_COMM_WORLD to initial method.
    bundle.put(cutoffRadiusFactor);
    bundle.put(latticeConst);

    // step and step length
    bundle.put(timeSteps);
    bundle.put(timeStepLength);

    bundle.put(createPhaseMode);
    bundle.put(createTSet);
    bundle.put(createSeed);
    bundle.put(readPhaseFilename);

    // alloy
    bundle.put(alloyCreateSeed);
    bundle.put(types.size());
    for (const AtomType &atom_ty:types) {
        atom_ty.packdata(bundle);
    }

    bundle.put(potentialFileType);
    bundle.put(potentialFilename);

    // output section
    bundle.put(output.presets.size());
    for (DumpConfig out: output.presets) {
        out.packdata(bundle);
    }

    bundle.put(output.thermo_interval);
    // logs subsection in output section.
    bundle.put(output.logs_mode);
    bundle.put(output.logs_filename);

    // stages
    bundle.put(stages.size());
    for (Stage stage:stages) {
        stage.packdata(bundle);
    }
}

void ConfigValues::unpackdata(kiwi::Bundle &bundle) {
    // fetch data from buffer.
//    if (getPackedData() != nullptr) { // buffer != null
    int cursor = 0;

    bundle.get(cursor, DIMENSION, phaseSpace);
    bundle.get(cursor, cutoffRadiusFactor);
    bundle.get(cursor, latticeConst);

    // step and step length
    bundle.get(cursor, timeSteps);
    bundle.get(cursor, timeStepLength);

    bundle.get(cursor, createPhaseMode);
    bundle.get(cursor, createTSet);
    bundle.get(cursor, createSeed);
    bundle.get(cursor, readPhaseFilename);

    // alloy
    bundle.get(cursor, alloyCreateSeed);
    std::size_t atom_types_size = 0;
    bundle.get(cursor, atom_types_size);
    types.resize(atom_types_size);
    for (std::size_t i = 0; i < atom_types_size; i++) {
        types[i].unnpackdata(cursor, bundle);
    }

    bundle.get(cursor, potentialFileType);
    bundle.get(cursor, potentialFilename);

    // output section.
    std::size_t presets_size = 0;
    bundle.get(cursor, presets_size);
    output.presets.resize(presets_size);
    for (std::size_t i = 0; i < presets_size; i++) {
        output.presets[i].unnpackdata(cursor, bundle);
    }

    bundle.get(cursor, output.thermo_interval);

    // logs subsection in output section.
    bundle.get(cursor, output.logs_mode);
    bundle.get(cursor, output.logs_filename);

    // stages
    std::size_t stages_size = 0;
    bundle.get(cursor, stages_size);
    stages.resize(stages_size);
    for (std::size_t i = 0; i < stages_size; i++) {
        stages[i].unnpackdata(cursor, bundle);
    }
}

std::ostream &operator<<(std::ostream &os, const ConfigValues &cv) {
    // simulation section
    os << "===========config of simulation=============" << std::endl;
    os << "simulation.phase_space:" << cv.phaseSpace[0] << "," << cv.phaseSpace[1]
       << "," << cv.phaseSpace[2] << "," << std::endl;
    os << "simulation.cutoff_radius:" << cv.cutoffRadiusFactor << std::endl;
    os << "simulation.lattice_const:" << cv.latticeConst << std::endl;
    os << "simulation.timesteps:" << cv.timeSteps << std::endl;

    os << "simulation.createphase.createPhaseMode:" << (cv.createPhaseMode ? "true" : "false") << std::endl;
    os << "simulation.createphase.tSet:" << cv.createTSet << std::endl;
    os << "simulation.createphase.seed:" << cv.createSeed << std::endl;
    os << "simulation.createphase.readPhaseFilename:" << cv.readPhaseFilename << std::endl;

    // simulation.alloy
    os << "simulation.alloy.alloyCreateSeed:" << cv.alloyCreateSeed << std::endl;
    os << "simulation.alloy.types:" << std::endl;
    for (const AtomType &atom_type : cv.types) {
        os << "Name: " << atom_type.name << ", mass: " << atom_type.mass << ", weight: " << atom_type.weight
           << std::endl;
    }
    os << "simulation.potential_file.type:" << cv.potentialFileType << std::endl;
    os << "simulation.potential_file.filename:" << cv.potentialFilename << std::endl;

    // output section
    os << "output.preset:" << std::endl;
    for (const DumpConfig &out : cv.output.presets) {
        os << "output.presets.name:" << out.name << std::endl;
        os << "output.presets.mode(debug:0, copy:1):" << out.mode << std::endl;
        os << "output.presets.region:" << out.region[0] << "," << out.region[1] << "," << out.region[2] << ","
           << out.region[3] << "," << out.region[4] << "," << out.region[5] << "," << std::endl;
        os << "output.presets.steps:" << out.steps << std::endl;
        os << "output.presets.file_path:" << out.file_path << std::endl;
        os << "output.presets.by-frame:" << out.by_frame << std::endl;
        os << "output.presets.dump-mask: " << out.dump_mask << std::endl;
        os << "output.presets.whole-system:" << out.dump_whole_system << std::endl;
    }
    os << "output.thermo.interval:" << cv.output.thermo_interval << std::endl;
    os << "output.logs.mode: "
       << (cv.output.logs_mode == LOGS_MODE_CONSOLE ? LOGS_MODE_CONSOLE_STRING : LOGS_MODE_FILE_STRING)
       << std::endl;
    os << "output.logs_filename:" << cv.output.logs_filename << std::endl;

    // stages
    os << "stages:" << std::endl;
    for (const Stage &stage : cv.stages) {
        os << "stage.steps:" << stage.steps << std::endl;
        os << "stage.steps_length:" << stage.step_length << std::endl;
        if (stage.collision_set) {
            os << "stage.collision.collision_step:" << stage.collisionStep << std::endl;
            os << "stage.collision.lat:" << stage.collisionLat[0] << "," << stage.collisionLat[1] << ","
               << stage.collisionLat[2] << "," << stage.collisionLat[3] << std::endl;
            os << "pka:" << stage.pkaEnergy << ",";
            os << "direction:" << stage.direction[0] << "," << stage.direction[1] << ","
               << stage.direction[2] << std::endl;
        }
        if (stage.rescales_set) {
            os << "rescale to T:" << stage.rescale_t << " every " << stage.rescale_every << " step(s)" << std::endl;
        }
        if (stage.dump_set) {
            os << "use dump preset:" << stage.dump_preset_use << " every " << stage.dump_every_steps << " step(s)"
               << std::endl;
        }
    }
    os << "============================================" << std::endl << std::endl;
    return os;
}
