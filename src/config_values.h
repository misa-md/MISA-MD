//
// Created by genshen(genshenchu@gmail.com) on 2018-3-9.
//

#ifndef CRYSTAL_MD_CONFIG_VALUES_H
#define CRYSTAL_MD_CONFIG_VALUES_H

#include <string>
#include <utils/bundle.h>
#include "pre_config.h"

using namespace std;

#define OUTPUT_COPY_MODE 0
#define OUTPUT_DIRECT_MODE 1
#define DEFAULT_OUTPUT_DUMP_FILENAME "crystal_md.out"

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
    string readPhaseFilename; // for read mode

    unsigned long collisionSteps;
    int collisionLat[4];
    double collisionV[DIMENSION];

    string potentialFileType;
    string potentialFilename;
    // simulation section ends
    // output section
    int outputMode;
    string outputDumpFilename;
    // output section ends
    // config values ends

    ConfigValues();

    void packdata(kiwi::Bundle &bundle);  // todo override

    void unpackdata(kiwi::Bundle &bundle);

};


#endif //CRYSTAL_MD_CONFIG_VALUES_H
