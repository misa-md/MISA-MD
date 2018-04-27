//
// Created by genshen(genshenchu@gmail.com) on 2017/4/16.
//
#include <toml.hpp>
#include "config_values.h"
#include "config/config.h"

#ifndef CRYSTAL_MD_CONFIG_PARSER_H
#define CRYSTAL_MD_CONFIG_PARSER_H

using namespace std;

/**
 * in our design, rank MASTER read config file,
 * and sync the config data to other processors.
 */
class ConfigParser : public kiwi::config {
public:
    ConfigParser();

    /**
     * Error flag while resolving toml config file.
     */
    bool hasError;
    /**
     * The error message  while resolving toml config file.
     */
    string errorMessage;

    // config values will be stored here.
    ConfigValues configValues;

    /**
     * Create a new config instance first,
     * then it will try to resolve the toml config file.
     * and config data from resolving will be stored in {@var configValues}
     * @param configureFilePath path of toml config file.
     * @return an pointer to {@class ConfigParser}.
     */
    static ConfigParser *newInstance(const string &configureFilePath);

    /**
     * Just create an empty {@class ConfigParser} class (not resolve toml config file).
     * @return
     */
    static ConfigParser *getInstance();

    bool configureCheck();

private:

    // pointer of this {@class ConfigParser} for single mode.
    static ConfigParser *m_pInstance; // stored in static area.

    /**
     * resolve toml format config using lib: https://github.com/skystrife/cpptoml
     * visit the url above to get more information.
     * @param table
     */
    void resolveConfig(std::shared_ptr<cpptoml::table> table) override;

    // resolve "simulation" section of toml config file.
    void resolveConfigSimulation(std::shared_ptr<cpptoml::table> v);

    // resolve "output" section of toml config file.
    void resolveConfigOutput(shared_ptr<cpptoml::table> v);

    /**
     * [master] put data into bundle, in which bundle is used to buffer config data.
     * later, the bundle will be sent to the other processors (sync config).
     * @param bundle used for buffering config data.
     */
    void putConfigData(kiwi::Bundle &bundle) override;

    /**
     * [non-master] get config data from bundle after finishing sync.
     * @param bundle
     */
    void getConfigData(kiwi::Bundle &bundle) override;
};


#endif //CRYSTAL_MD_CONFIG_PARSER_H
