//
// Created by genshen on 2018-3-9.
//

#ifndef CRYSTAL_MD_CONFIG_VALUES_H
#define CRYSTAL_MD_CONFIG_VALUES_H

#include <string>
#include "utils/data_pack.h"
#include "pre_config.h"

using namespace std;

class ConfigValues : public DataPack {
    friend std::ostream &operator<<(std::ostream &os, const ConfigValues &cv);

public:
    // simulation section
    int phaseSpace[DIMENSION];
    double cutoffRadius;
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

    ConfigValues(unsigned int cap);

//    ~ConfigValues(); //todo remove
    void packdata();  // todo override

    void unpackdata();

};


#endif //CRYSTAL_MD_CONFIG_VALUES_H
