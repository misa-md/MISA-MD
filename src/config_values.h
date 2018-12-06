//
// Created by genshen(genshenchu@gmail.com) on 2018-3-9.
//

#ifndef CRYSTAL_MD_CONFIG_VALUES_H
#define CRYSTAL_MD_CONFIG_VALUES_H

#include <string>
#include <utils/bundle.h>
#include "types/pre_define.h"
#include "types/atom_types.h"

#define OUTPUT_COPY_MODE 0
#define OUTPUT_DIRECT_MODE 1
#define LOGS_MODE_CONSOLE 0
#define LOGS_MODE_FILE 1
#define LOGS_MODE_CONSOLE_STRING "console"
#define LOGS_MODE_FILE_STRING "file"
#define DEFAULT_LOGS_MODE_CONSOLE_STRING LOGS_MODE_CONSOLE_STRING
#define DEFAULT_OUTPUT_DUMP_FILE_PATH "crystal_md.out"

typedef short _type_out_mode;
typedef short _type_logs_mode;

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
    double collisionV[DIMENSION];

    std::string potentialFileType;
    std::string potentialFilename;
    // simulation section ends
    // output section
    _type_out_mode atomsDumpMode;
    unsigned long atomsDumpInterval;
    bool outByFrame; // output atoms by frame if true.
    std::string atomsDumpFilePath;
    // logs in output section
    _type_logs_mode logs_mode;
    std::string logs_filename;
    // output section ends
    // config values ends

    ConfigValues();

    void packdata(kiwi::Bundle &bundle);  // todo override

    void unpackdata(kiwi::Bundle &bundle);

};


#endif //CRYSTAL_MD_CONFIG_VALUES_H
