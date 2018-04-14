//
// Created by genshen(genshenchu@gmail.com) on 2017/4/16.
//
#include <toml.hpp>
#include "config_values.h"
#include "config/config.h"

#ifndef CRYSTAL_MD_CONFIG_H
#define CRYSTAL_MD_CONFIG_H

using namespace std;

class Config : public kiwi::config {
public:
    Config();

    bool hasError;
    string errorMessage;

    // config values
    ConfigValues configValues;

    static Config *newInstance(const string &configureFilePath);

    static Config *getInstance();

    bool configureCheck();

private:

    // when sync, some data like string is first copied here; then use MPI_Bcast to send whole class to other processors,
    // and the other processors resume data (e.g. string) from buffer array.
//    char buffer[512];

    static Config *m_pInstance; // stored in static area.

//    config(const string &configurePath);

    void resolveConfig(std::shared_ptr<cpptoml::table> table) override;

    // resolve "simulation" section of config file.
    void resolveConfigSimulation(std::shared_ptr<cpptoml::table> v);

    // resolve "output" section of config file.
    void resolveConfigOutput(shared_ptr<cpptoml::table> v);

    void putConfigData(kiwi::Bundle &bundle) override;

    void getConfigData(kiwi::Bundle &bundle) override;
};


#endif //CRYSTAL_MD_CONFIG_H
