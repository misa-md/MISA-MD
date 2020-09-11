//
// Created by genshen(genshenchu@gmail.com) on 2018-3-9.
//

#ifndef MISA_MD_CONFIG_VALUES_H
#define MISA_MD_CONFIG_VALUES_H

#include <string>
#include <vector>
#include <utils/bundle.h>
#include "types/pre_define.h"
#include "types/atom_types.h"

enum OutputMode {
    DEBUG = 0,
    COPY = 1,
};

#define LOGS_MODE_CONSOLE 0
#define LOGS_MODE_FILE 1
#define LOGS_MODE_CONSOLE_STRING "console"
#define LOGS_MODE_FILE_STRING "file"
#define DEFAULT_LOGS_MODE_CONSOLE_STRING LOGS_MODE_CONSOLE_STRING
#define DEFAULT_OUTPUT_DUMP_FILE_PATH "crystal_md.out"
#define ORIGIN_OUTPUT_DUMP_FILE_PATH "origin_crystal_md.out"

typedef short _type_logs_mode;

struct Stage {
    unsigned long steps;
    double step_length;

    // collision
    bool collision_set;
    unsigned long collisionStep;
    int collisionLat[4];
    double pkaEnergy;
    double direction[DIMENSION];

    // rescale
    bool rescales_set;
    double rescale_t; // rescale to a temperature.
    unsigned long rescale_every; // step to do rescale

    Stage();

    void packdata(kiwi::Bundle &bundle);

    void unnpackdata(int &cursor, kiwi::Bundle &bundle);
};

struct Output {
    // output section
    OutputMode atomsDumpMode;
    unsigned long atomsDumpInterval;
    // output atoms by frame if true.
    bool outByFrame;
    std::string atomsDumpFilePath;
    // path of dumped origin atoms before collision
    std::string originDumpPath;
    // interval to output thermodynamics information.
    unsigned long thermo_interval;
    // logs in output section
    _type_logs_mode logs_mode;
    std::string logs_filename;

    Output() : atomsDumpMode(OutputMode::COPY), atomsDumpInterval(1),
               outByFrame(false), originDumpPath(ORIGIN_OUTPUT_DUMP_FILE_PATH),
               atomsDumpFilePath(DEFAULT_OUTPUT_DUMP_FILE_PATH), thermo_interval(0),
               logs_mode(LOGS_MODE_CONSOLE), logs_filename("") {}
};

class ConfigValues {
    friend std::ostream &operator<<(std::ostream &os, const ConfigValues &cv);

public:
    // config values start
    // simulation section
    int64_t phaseSpace[DIMENSION];
    double cutoffRadiusFactor;
    double latticeConst;
    unsigned long timeSteps; // total steps is not set in config file, but compute from each stages.
    double timeStepLength; // default step length

    bool createPhaseMode;
    double createTSet; // system temperature for creation.
    int createSeed;
    std::string readPhaseFilename; // for read mode

    // alloy
    int alloyCreateSeed;
    int alloyRatio[atom_type::num_atom_types];

    // potential config
    std::string potentialFileType;
    std::string potentialFilename;
    // simulation section ends
    // output section
    Output output;
    // config values ends
    std::vector<Stage> stages;

    ConfigValues();

    void packdata(kiwi::Bundle &bundle);  // todo override

    void unpackdata(kiwi::Bundle &bundle);

};


#endif //MISA_MD_CONFIG_VALUES_H
