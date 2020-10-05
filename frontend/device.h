//
// Created by genshen on 2019-05-30.
//

#ifndef MISA_MD_DEVICE_H
#define MISA_MD_DEVICE_H

#include <unistd.h>   // for isatty()
#include <stdio.h>    // for fileno()

#ifdef __cplusplus
extern "C"
{
#endif

int istty() {
    return isatty(fileno(stdout));
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif //MISA_MD_DEVICE_H
