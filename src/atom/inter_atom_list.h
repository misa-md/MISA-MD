//
// Created by genshen on 2018-05-19.
//

#ifndef CRYSTAL_MD_INTER_ATOM_LIST_H
#define CRYSTAL_MD_INTER_ATOM_LIST_H

#include <vector>
#include "../pre_define.h"
#include "atom_types.h"

/**
 * storing inter atoms
 */
class InterAtomList {
public:
    std::vector<_type_atom_id > idinter;
    std::vector<atom_type::atom_type > typeinter;
    std::vector<std::vector<double>> xinter; // 间隙原子坐标
    std::vector<std::vector<double>> vinter; // 间隙原子速度
    std::vector<std::vector<double>> finter; // 间隙原子力
    std::vector<double> rhointer;
    std::vector<double> dfinter;
    int nlocalinter, nghostinter; // 本地间隙原子数和ghost间隙原子数

    InterAtomList();

    void appendInter(_type_atom_id atom_id);
};


#endif //CRYSTAL_MD_INTER_ATOM_LIST_H
