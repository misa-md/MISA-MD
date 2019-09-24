//
// Created by genshen(genshenchu@gmail.com) on 2018-3-9.
//

#ifndef CRYSTAL_MD_CONFIG_VALUES_H
#define CRYSTAL_MD_CONFIG_VALUES_H

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
    unsigned long timeSteps;

    double timeStepLength;
    size_t vsl_size = 0; // array size of variable time step length.
    std::vector<unsigned long> vsl_break_points;
    std::vector<double> vsl_lengths;

    bool createPhaseMode;
    double createTSet;
    int createSeed;
    std::string readPhaseFilename; // for read mode

    // alloy
    int alloyCreateSeed;
    int alloyRatio[atom_type::num_atom_types];

    // collision
    unsigned long collisionStep;
    int collisionLat[4];
    double pkaEnergy;
    double direction[DIMENSION];

    std::string potentialFileType;
    std::string potentialFilename;
    // simulation section ends
    // output section
    Output output;
    // config values ends

    ConfigValues();

    void setVarStepLengths(const unsigned long *break_points, const double *lengths, const unsigned long size);

    /**
     * set variable time step length config.
     * @param break_points break points of step
     * @param lengths time step length after each break pointer step.
     * @param size size of array @param break_points and @param lengths.
     */
    void setVarStepLengths(const std::vector<unsigned long> break_points, const std::vector<double> lengths,
                           const unsigned long size);

    void packdata(kiwi::Bundle &bundle);  // todo override

    void unpackdata(kiwi::Bundle &bundle);

};


#endif //CRYSTAL_MD_CONFIG_VALUES_H
