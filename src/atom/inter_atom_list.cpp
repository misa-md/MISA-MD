//
// Created by genshen on 5/19/18.
//

#include "inter_atom_list.h"

InterAtomList::InterAtomList() : nlocalinter(0), nghostinter(0) {

}

void InterAtomList::appendInter(_type_atom_id atom_id) {

}

void InterAtomList::addInterAtom(AtomElement &atom) {
    inter_list.push_back(atom);
    nlocalinter++;
    // todo set df,f,rhointer to 0.
}

void InterAtomList::pack_intersend(std::vector<unsigned long> interbuf, particledata *buf) {
    int j;
    for (int i = 0; i < interbuf.size(); i++) {
        j = interbuf[i];
        buf[i].id = inter->idinter[j];
        buf[i].type = inter->typeinter[j];
        buf[i].r[0] = inter->xinter[j][0];
        buf[i].r[1] = inter->xinter[j][1];
        buf[i].r[2] = inter->xinter[j][2];
        buf[i].v[0] = inter->vinter[j][0];
        buf[i].v[1] = inter->vinter[j][1];
        buf[i].v[2] = inter->vinter[j][2];
        // remove the inter atom.
        // exchange the atom at end of vector to the position of atom (inter[j]) te be removed .
        inter->idinter[j] = inter->idinter[inter->nlocalinter - 1];
        inter->typeinter[j] = inter->typeinter[inter->nlocalinter - 1];
        inter->xinter[j][0] = inter->xinter[inter->nlocalinter - 1][0];
        inter->xinter[j][1] = inter->xinter[inter->nlocalinter - 1][1];
        inter->xinter[j][2] = inter->xinter[inter->nlocalinter - 1][2];
        inter->vinter[j][0] = inter->vinter[inter->nlocalinter - 1][0];
        inter->vinter[j][1] = inter->vinter[inter->nlocalinter - 1][1];
        inter->vinter[j][2] = inter->vinter[inter->nlocalinter - 1][2];
        inter->nlocalinter--;
    }
}

void InterAtomList::unpack_interrecv(int d, int n,
                                     double lower[DIMENSION], // p_domain->getMeasuredSubBoxLowerBounding(d)
                                     double upper[DIMENSION], // p_domain->getMeasuredSubBoxUpperBounding(d)
                                     particledata *buf) {
    std::vector<double> xtemp(3);
    std::vector<double> vtemp(3);
    unsigned long id;
    atom_type::atom_type type;
    for (int i = 0; i < n; i++) {
        id = buf[i].id;
        type = buf[i].type;
        xtemp[0] = buf[i].r[0];
        xtemp[1] = buf[i].r[1];
        xtemp[2] = buf[i].r[2];
        vtemp[0] = buf[i].v[0];
        vtemp[1] = buf[i].v[1];
        vtemp[2] = buf[i].v[2];
        if (xtemp[d] >= lower[d] &&
            xtemp[d] < upper[d]) {
            if (inter->nlocalinter == inter->xinter.size()) {
                inter->idinter.push_back(id);
                inter->typeinter.push_back(type);
                inter->xinter.push_back(xtemp);
                inter->vinter.push_back(vtemp);
                inter->nlocalinter++;
                inter->finter.resize(inter->nlocalinter, std::vector<double>(3));
                inter->rhointer.resize(inter->nlocalinter);
                inter->dfinter.resize(inter->nlocalinter);
            } else {
                if (inter->idinter.size() == inter->nlocalinter) {
                    inter->idinter.push_back(id);
                } else {
                    inter->idinter[inter->nlocalinter] = id;
                }
                inter->typeinter[inter->nlocalinter] = type;
                inter->xinter[inter->nlocalinter][0] = xtemp[0];
                inter->xinter[inter->nlocalinter][1] = xtemp[1];
                inter->xinter[inter->nlocalinter][2] = xtemp[2];
                if (inter->nlocalinter == inter->vinter.size()) {
                    inter->vinter.push_back(vtemp);
                } else {
                    inter->vinter[inter->nlocalinter][0] = vtemp[0];
                    inter->vinter[inter->nlocalinter][1] = vtemp[1];
                    inter->vinter[inter->nlocalinter][2] = vtemp[2];
                }
                inter->nlocalinter++;
                inter->finter.resize(inter->nlocalinter, std::vector<double>(3));
                inter->rhointer.resize(inter->nlocalinter);
                inter->dfinter.resize(inter->nlocalinter);
            }
        }
    }
}

void InterAtomList::pack_bordersend(int dimension, int n,
                                    std::vector<int> &sendlist, LatParticleData *buf, double shift) {
    int j;
    if (dimension == 0) {
        for (int i = 0; i < n; i++) {
            j = sendlist[i];
            buf[i].type = inter->typeinter[j];
            buf[i].r[0] = inter->xinter[j][0] + shift;
            buf[i].r[1] = inter->xinter[j][1];
            buf[i].r[2] = inter->xinter[j][2];
        }
    } else if (dimension == 1) {
        for (int i = 0; i < n; i++) {
            j = sendlist[i];
            buf[i].type = inter->typeinter[j];
            buf[i].r[0] = inter->xinter[j][0];
            buf[i].r[1] = inter->xinter[j][1] + shift;
            buf[i].r[2] = inter->xinter[j][2];
        }
    } else {
        for (int i = 0; i < n; i++) {
            j = sendlist[i];
            buf[i].type = inter->typeinter[j];
            buf[i].r[0] = inter->xinter[j][0];
            buf[i].r[1] = inter->xinter[j][1];
            buf[i].r[2] = inter->xinter[j][2] + shift;
        }
    }
}

void InterAtomList::unpack_borderrecv(int n,
                                      double lower[DIMENSION], // p_domain->getMeasuredGhostLowerBounding(d)
                                      double upper[DIMENSION], // p_domain->getMeasuredGhostUpperBounding(d)
                                      LatParticleData *buf, std::vector<int> &recvlist) {
    atom_type::atom_type type;
    std::vector<double> xtemp(3);
    for (int i = 0; i < n; i++) {
        type = buf[i].type;
        xtemp[0] = buf[i].r[0];
        xtemp[1] = buf[i].r[1];
        xtemp[2] = buf[i].r[2];
        if (xtemp[0] >= lower[0] &&
            xtemp[0] < upper[0] &&
            xtemp[1] >= lower[1] &&
            xtemp[1] < upper[1] &&
            xtemp[2] >= lower[2] &&
            xtemp[2] < upper[2]) {
            if (inter->xinter.size() == inter->nlocalinter + inter->nghostinter) {
                inter->typeinter.push_back(type);
                inter->xinter.push_back(xtemp);
                inter->nghostinter++;
                recvlist[i] = inter->nlocalinter + inter->nghostinter - 1;
                inter->finter.resize(inter->nlocalinter + inter->nghostinter, std::vector<double>(3));
                inter->rhointer.resize(inter->nlocalinter + inter->nghostinter);
                inter->dfinter.resize(inter->nlocalinter + inter->nghostinter);
            } else {
                inter->typeinter[inter->nlocalinter + inter->nghostinter] = type;
                inter->xinter[inter->nlocalinter + inter->nghostinter][0] = xtemp[0];
                inter->xinter[inter->nlocalinter + inter->nghostinter][1] = xtemp[1];
                inter->xinter[inter->nlocalinter + inter->nghostinter][2] = xtemp[2];
                inter->nghostinter++;
                recvlist[i] = inter->nlocalinter + inter->nghostinter - 1;
                inter->finter.resize(inter->nlocalinter + inter->nghostinter, std::vector<double>(3));
                inter->rhointer.resize(inter->nlocalinter + inter->nghostinter);
                inter->dfinter.resize(inter->nlocalinter + inter->nghostinter);
            }
        } else {
            recvlist[i] = -1;
        }
    }
}
