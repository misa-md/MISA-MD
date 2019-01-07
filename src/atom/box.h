//
// Created by genshen on 2019-01-07.
//

#ifndef CRYSTAL_MD_BOX_H
#define CRYSTAL_MD_BOX_H

/**
 * box status of atoms
 */
namespace box {
    typedef unsigned int _type_flag_32;
    const _type_flag_32 IN_BOX = 0; // in box
    const _type_flag_32 OUT_BOX_X_LITTER = 1; // out of box at x direction(in litter end)
    const _type_flag_32 OUT_BOX_X_BIG = 1 << 1; // out of box at x direction(in big end)
    const _type_flag_32 OUT_BOX_Y_LITTER = 1 << 2; // out of box at y direction(in litter end)
    const _type_flag_32 OUT_BOX_Y_BIG = 1 << 3; // out of box at y direction(in big end)
    const _type_flag_32 OUT_BOX_Z_LITTER = 1 << 4; // out of box at z direction(in litter end)
    const _type_flag_32 OUT_BOX_Z_BIG = 1 << 5; // out of box at z direction(in big end)
}

#endif //CRYSTAL_MD_BOX_H
