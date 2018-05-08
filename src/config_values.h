//
// Created by genshen(genshenchu@gmail.com) on 2018-3-9.
//

#ifndef CRYSTAL_MD_CONFIG_VALUES_H
#define CRYSTAL_MD_CONFIG_VALUES_H

#include <string>
#include <utils/bundle.h>
#include "pre_define.h"
#include "atom_types.h"

#define OUTPUT_COPY_MODE 0
#define OUTPUT_DIRECT_MODE 1
#define DEFAULT_OUTPUT_DUMP_FILENAME "crystal_md.out"

typedef short _type_out_mode;

class ConfigValues {
    friend std::ostream &operator<<(std::ostream &os, const ConfigValues &cv);

public:
    // config values start
    // simulation section
    int64_t phaseSpace[DIMENSION];
    double cutoffRadiusFactor;
    double latticeConst;
    unsigned long timeSteps;

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
    _type_out_mode outputMode;
    std::string outputDumpFilename;
    // output section ends
    // config values ends

    ConfigValues();

    void packdata(kiwi::Bundle &bundle);  // todo override

    void unpackdata(kiwi::Bundle &bundle);

};


#endif //CRYSTAL_MD_CONFIG_VALUES_H
