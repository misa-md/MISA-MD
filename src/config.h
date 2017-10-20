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

    bool hasError = false;
    string errorMessage;

    //config values
    int phaseSpace[DIMENSION] = {0, 0, 0};
    double cutoffRadius = 0.0;
    double latticeConst = 0.0;
    unsigned long timeSteps = 100;

    bool createPhaseMode = true;
    double createTSet = 0.0;
    int createSeed = 0;
    string readPhaseFilename = ""; // for read mode

    unsigned long collisionSteps = 0;
    int collisionLat[4] = {0, 0, 0, 0};
    double collisionV[DIMENSION] = {0.0, 0.0, 0.0};

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
