//
// Created by genshen(genshenchu@gmail.com) on 2017/4/16.
//
#include "toml.hpp"
#include "config_values.h"

#ifndef CRYSTAL_MD_CONFIG_H
#define CRYSTAL_MD_CONFIG_H

using namespace std;

class config {
public:
    config();

    bool hasError;
    string errorMessage;

    // config values
    ConfigValues *configValues;

    static config *newInstance(const string &configureFilePath);

    static config *newInstance();

    bool configureCheck();

    void sync(); // Synchronize configure information to all processors.

private:

    // when sync, some data like string is first copied here; then use MPI_Bcast to send whole class to other processors,
    // and the other processors resume data (e.g. string) from buffer array.
//    char buffer[512];

    static config *m_pInstance; // stored in static area.

    config(const string &configurePath);

    void resolveConfig(const string &configurePath);

    // resolve "simulation" section of config file.
    void resolveConfigSimulation(const toml::Value &v);

    // resolve "output" section of config file.
    void resolveConfigOutput(const toml::Value &v);
};


#endif //CRYSTAL_MD_CONFIG_H
