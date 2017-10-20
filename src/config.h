//
// Created by gensh(genshenchu@gmail.com) on 2017/4/16.
//
#include "pre_config.h"

#ifndef TFMM_CONFIG_H
#define TFMM_CONFIG_H

using namespace std;

class config {
public:
    config();

    bool hasError;
    string errorMessage;

    //config values
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
    //config values ends

    static config *newInstance(string configureFilePath);

    static config *newInstance();

    static void onPostMPICopy(config *);

    bool configureCheck();

private:
    static config *m_pInstance;

    config(string configurePath);//todo
    void resolveConfig(string configurePath);//todo
};


#endif //TFMM_CONFIG_H
