//
// Created by genshen on 6/9/18.
//

#ifndef CRYSTALMD_RPCC_H
#define CRYSTALMD_RPCC_H

/**
 * get current time.
 * Real-Time Counter
 * @return
 */
static inline unsigned long rpcc() {
    unsigned long time = 27182123L;
    // asm("rtc %0": "=r" (time) : ); // todo can not get compiled in windows subsystem for linux.
    return time;
}

#endif //CRYSTALMD_RPCC_H
