//
// Created by genshen on 5/17/18.
//

#include <logs/logs.h>
#include "pack.h"
#include "../atom_element.h"

/*
_atom->nlocalinter,
//_atom->nghostinter,
_atom->idinter,
_atom->typeinter,
_atom->xinter,
_atom->vinter,
_atom->finter,
_atom->rhointer,
_atom->dfinter,
 */

void pack::pack_intersend(
        int &nlocalinter,
//        int &nghostinter,
        std::vector<_type_atom_id> &idinter,
        std::vector<_type_atom_type> &typeinter,
        std::vector<std::vector<double>> &xinter,
        std::vector<std::vector<double>> &vinter,
//        std::vector<std::vector<double>> &finter,
//        std::vector<double> &rhointer,
//        std::vector<double> &dfinter,
        std::vector<unsigned long> interbuf, particledata *buf) {
    int j;
    for (int i = 0; i < interbuf.size(); i++) {
        j = interbuf[i];
        buf[i].id = idinter[j];
        buf[i].type = typeinter[j];
        buf[i].r[0] = xinter[j][0];
        buf[i].r[1] = xinter[j][1];
        buf[i].r[2] = xinter[j][2];
        buf[i].v[0] = vinter[j][0];
        buf[i].v[1] = vinter[j][1];
        buf[i].v[2] = vinter[j][2];
        // remove the inter atom.
        // exchange the atom at end of vector to the position of atom (inter[j]) te be removed .
        idinter[j] = idinter[nlocalinter - 1];
        typeinter[j] = typeinter[nlocalinter - 1];
        xinter[j][0] = xinter[nlocalinter - 1][0];
        xinter[j][1] = xinter[nlocalinter - 1][1];
        xinter[j][2] = xinter[nlocalinter - 1][2];
        vinter[j][0] = vinter[nlocalinter - 1][0];
        vinter[j][1] = vinter[nlocalinter - 1][1];
        vinter[j][2] = vinter[nlocalinter - 1][2];
        nlocalinter--;
    }
}

void pack::unpack_interrecv(int d, int n, int &nlocalinter,
                            std::vector<_type_atom_id> &idinter,
                            std::vector<_type_atom_type> &typeinter,
                            std::vector<std::vector<double>> &xinter,
                            std::vector<std::vector<double>> &vinter,
                            std::vector<std::vector<double>> &finter,
                            std::vector<double> &rhointer,
                            std::vector<double> &dfinter,
                            double lower[DIMENSION], // p_domain->getMeasuredSubBoxLowerBounding(d)
                            double upper[DIMENSION], // p_domain->getMeasuredSubBoxUpperBounding(d)
                            particledata *buf) {
    std::vector<double> xtemp(3);
    std::vector<double> vtemp(3);
    unsigned long id;
    int type;
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
            if (nlocalinter == xinter.size()) {
                idinter.push_back(id);
                typeinter.push_back(type);
                xinter.push_back(xtemp);
                vinter.push_back(vtemp);
                nlocalinter++;
                finter.resize(nlocalinter, std::vector<double>(3));
                rhointer.resize(nlocalinter);
                dfinter.resize(nlocalinter);
            } else {
                if (idinter.size() == nlocalinter) {
                    idinter.push_back(id);
                } else {
                    idinter[nlocalinter] = id;
                }
                typeinter[nlocalinter] = type;
                xinter[nlocalinter][0] = xtemp[0];
                xinter[nlocalinter][1] = xtemp[1];
                xinter[nlocalinter][2] = xtemp[2];
                if (nlocalinter == vinter.size()) {
                    vinter.push_back(vtemp);
                } else {
                    vinter[nlocalinter][0] = vtemp[0];
                    vinter[nlocalinter][1] = vtemp[1];
                    vinter[nlocalinter][2] = vtemp[2];
                }
                nlocalinter++;
                finter.resize(nlocalinter, std::vector<double>(3));
                rhointer.resize(nlocalinter);
                dfinter.resize(nlocalinter);
            }
        }
    }
}

void pack::pack_bordersend(int dimension, int n,
                           std::vector<_type_atom_type> &typeinter,
                           std::vector<std::vector<double>> &xinter,
                           std::vector<int> &sendlist, LatParticleData *buf, double shift) {
    int j;
    if (dimension == 0) {
        for (int i = 0; i < n; i++) {
            j = sendlist[i];
            buf[i].type = typeinter[j];
            buf[i].r[0] = xinter[j][0] + shift;
            buf[i].r[1] = xinter[j][1];
            buf[i].r[2] = xinter[j][2];
        }
    } else if (dimension == 1) {
        for (int i = 0; i < n; i++) {
            j = sendlist[i];
            buf[i].type = typeinter[j];
            buf[i].r[0] = xinter[j][0];
            buf[i].r[1] = xinter[j][1] + shift;
            buf[i].r[2] = xinter[j][2];
        }
    } else {
        for (int i = 0; i < n; i++) {
            j = sendlist[i];
            buf[i].type = typeinter[j];
            buf[i].r[0] = xinter[j][0];
            buf[i].r[1] = xinter[j][1];
            buf[i].r[2] = xinter[j][2] + shift;
        }
    }
}

void pack::unpack_borderrecv(int n,
                             int &nlocalinter,
                             int &nghostinter,
//                             std::vector<_type_atom_id> &idinter,
                             std::vector<_type_atom_type> &typeinter,
                             std::vector<std::vector<double>> &xinter,
//                             std::vector<std::vector<double>> &vinter,
                             std::vector<std::vector<double>> &finter,
                             std::vector<double> &rhointer,
                             std::vector<double> &dfinter,
                             double lower[DIMENSION], // p_domain->getMeasuredGhostLowerBounding(d)
                             double upper[DIMENSION], // p_domain->getMeasuredGhostUpperBounding(d)
                             LatParticleData *buf, std::vector<int> &recvlist) {
    int type;
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
            if (xinter.size() == nlocalinter + nghostinter) {
                typeinter.push_back(type);
                xinter.push_back(xtemp);
                nghostinter++;
                recvlist[i] = nlocalinter + nghostinter - 1;
                finter.resize(nlocalinter + nghostinter, std::vector<double>(3));
                rhointer.resize(nlocalinter + nghostinter);
                dfinter.resize(nlocalinter + nghostinter);
            } else {
                typeinter[nlocalinter + nghostinter] = type;
                xinter[nlocalinter + nghostinter][0] = xtemp[0];
                xinter[nlocalinter + nghostinter][1] = xtemp[1];
                xinter[nlocalinter + nghostinter][2] = xtemp[2];
                nghostinter++;
                recvlist[i] = nlocalinter + nghostinter - 1;
                finter.resize(nlocalinter + nghostinter, std::vector<double>(3));
                rhointer.resize(nlocalinter + nghostinter);
                dfinter.resize(nlocalinter + nghostinter);
            }
        } else
            recvlist[i] = -1;
    }
}

/**
 *
 * In @var shift[i] ,i =0,1,2; at most one dimension has a non-zero shift[i].
 */
void pack::pack_send(const int dimension, const int n, const double shift[DIMENSION], AtomList &atom_list,
                     LatParticleData *buf, std::vector<_type_atom_id> &sendlist) {
    _type_atom_id j;
    for (int i = 0; i < n; i++) {
        j = sendlist[i];
        AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
        // for ghost atoms, we just care their position and atom type(eam calculating), so positions and types are enough.
        buf[i].type = atom_.type; // fixme --type
        buf[i].r[0] = atom_.x[0] + shift[0];
        buf[i].r[1] = atom_.x[1] + shift[1];
        buf[i].r[2] = atom_.x[2] + shift[2];
    }
}

void pack::unpack_recvfirst(int d, int direction, int n, AtomList &atom_list,
                            _type_lattice_size ghost[DIMENSION], //p_domain->getGhostLatticeSize(d)
                            _type_lattice_size box[DIMENSION], //p_domain->getSubBoxLatticeSize(d)
                            _type_lattice_size ext[DIMENSION], //p_domain->getGhostExtLatticeSize(d)
                            LatParticleData *buf, std::vector<std::vector<_type_atom_id> > &recvlist) {
    int xstart, ystart, zstart;
    int xstop, ystop, zstop;
    int recv_idnex = 2 * d + direction;
    int m = 0;
    if (d == 0) {
        if (direction == 0) { // mirror with send.
            xstart = ghost[0] + box[0];
            xstop = ext[0];
            ystart = ghost[1];
            ystop = ystart + box[1];
            zstart = ghost[2];
            zstop = zstart + box[2];
            for (int k = zstart; k < zstop; k++) {
                for (int j = ystart; j < ystop; j++) {
                    for (int i = xstart; i < xstop; i++) {
//                        kk = IndexOf3DIndex(i, j, k);
                        AtomElement &atom_ = atom_list.getAtomEleByGhostIndex(i, j, k);
                        atom_.type = buf[m].type;
                        atom_.x[0] = buf[m].r[0];
                        atom_.x[1] = buf[m].r[1];
                        atom_.x[2] = buf[m++].r[2];
                        recvlist[0].push_back(atom_list.IndexOf3DIndex(i,j,k));
                    }
                }
            }
            if (n != recvlist[0].size()) { // todo error.
                kiwi::logs::e("unpack_recvfirst", "received data size does not match the Mpi_Proble size.\n");
            }
        } else {
            xstart = 0;
            xstop = ghost[0];
            ystart = ghost[1];
            ystop = ystart + box[1];
            zstart = ghost[2];
            zstop = zstart + box[2];
            for (int k = zstart; k < zstop; k++) {
                for (int j = ystart; j < ystop; j++) {
                    for (int i = xstart; i < xstop; i++) {
//                        kk = IndexOf3DIndex(i, j, k);
                        AtomElement &atom_ = atom_list.getAtomEleByGhostIndex(i, j, k);
                        atom_.type = buf[m].type;
                        atom_.x[0] = buf[m].r[0];
                        atom_.x[1] = buf[m].r[1];
                        atom_.x[2] = buf[m++].r[2];
                        recvlist[1].push_back(atom_list.IndexOf3DIndex(i,j,k));
                    }
                }
            }
            if (n != recvlist[1].size()) // todo error handling in dataReuse feature.
                printf("wrong-01 !!!\n");
        }
    } else if (d == 1) {
        if (direction == 0) {
            xstart = 0;
            xstop = ext[0];
            ystart = ghost[1] + box[1];
            ystop = ext[1];
            zstart = ghost[2];
            zstop = zstart + box[2];
            for (int k = zstart; k < zstop; k++) {
                for (int j = ystart; j < ystop; j++) {
                    for (int i = xstart; i < xstop; i++) {
//                        kk = IndexOf3DIndex(i, j, k);
                        AtomElement &atom_ = atom_list.getAtomEleByGhostIndex(i, j, k);
                        atom_.type = buf[m].type;
                        atom_.x[0] = buf[m].r[0];
                        atom_.x[1] = buf[m].r[1];
                        atom_.x[2] = buf[m++].r[2];
                        recvlist[2].push_back(atom_list.IndexOf3DIndex(i,j,k));
                    }
                }
            }
            if (n != recvlist[2].size()) // todo error handling in dataReuse feature.
                printf("wrong-10!!!\n");
        } else {
            xstart = 0;
            xstop = ext[0];
            ystart = 0;
            ystop = ghost[1];
            zstart = ghost[2];
            zstop = zstart + box[2];
            for (int k = zstart; k < zstop; k++) {
                for (int j = ystart; j < ystop; j++) {
                    for (int i = xstart; i < xstop; i++) {
                        AtomElement &atom_ = atom_list.getAtomEleByGhostIndex(i, j, k);
                        atom_.type = buf[m].type;
                        atom_.x[0] = buf[m].r[0];
                        atom_.x[1] = buf[m].r[1];
                        atom_.x[2] = buf[m++].r[2];
                        recvlist[3].push_back(atom_list.IndexOf3DIndex(i,j,k));
                    }
                }
            }
            if (n != recvlist[3].size()) // todo error handling in dataReuse feature.
                printf("wrong-11 !!!\n");
        }
    } else {
        if (direction == 0) {
            xstart = 0;
            xstop = ext[0];
            ystart = 0;
            ystop = ext[1];
            zstart = ghost[2] + box[2];
            zstop = ext[2];
            for (int k = zstart; k < zstop; k++) {
                for (int j = ystart; j < ystop; j++) {
                    for (int i = xstart; i < xstop; i++) {
//                        kk = IndexOf3DIndex(i, j, k);
                        AtomElement &atom_ = atom_list.getAtomEleByGhostIndex(i, j, k);
                        atom_.type = buf[m].type;
                        atom_.x[0] = buf[m].r[0];
                        atom_.x[1] = buf[m].r[1];
                        atom_.x[2] = buf[m++].r[2];
                        recvlist[4].push_back(atom_list.IndexOf3DIndex(i,j,k));
                    }
                }
            }
            if (n != recvlist[4].size()) { // todo error handling in dataReuse feature.
                printf("wrong-20 !!\n");
            }
        } else {
            xstart = 0;
            xstop = ext[0];
            ystart = 0;
            ystop = ext[1];
            zstart = 0;
            zstop = ghost[2];
            for (int k = zstart; k < zstop; k++) {
                for (int j = ystart; j < ystop; j++) {
                    for (int i = xstart; i < xstop; i++) {
//                        kk = IndexOf3DIndex(i, j, k);
                        AtomElement &atom_ = atom_list.getAtomEleByGhostIndex(i, j, k);
                        atom_.type = buf[m].type;
                        atom_.x[0] = buf[m].r[0];
                        atom_.x[1] = buf[m].r[1];
                        atom_.x[2] = buf[m++].r[2];
                        recvlist[5].push_back(atom_list.IndexOf3DIndex(i,j,k));
                    }
                }
            }
            if (n != recvlist[5].size()) { // todo error handling in dataReuse feature.
                printf("wrong-30 !!!\n");
            }
        }
    }
}

void pack::unpack_recv(int d, int direction, int n, AtomList &atom_list,
                       LatParticleData *buf, std::vector<std::vector<_type_atom_id> > &recvlist) {
    long kk;
    if (d == 0) {
        if (direction == 0) {
            for (int i = 0; i < n; i++) {
                kk = recvlist[0][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(kk);
                atom_.type = buf[i].type;
                atom_.x[0] = buf[i].r[0];
                atom_.x[1] = buf[i].r[1];
                atom_.x[2] = buf[i].r[2];
            }
        } else {
            for (int i = 0; i < n; i++) {
                kk = recvlist[1][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(kk);
                atom_.type = buf[i].type;
                atom_.x[0] = buf[i].r[0];
                atom_.x[1] = buf[i].r[1];
                atom_.x[2] = buf[i].r[2];
            }
        }
    } else if (d == 1) {
        if (direction == 0) {
            for (int i = 0; i < n; i++) {
                kk = recvlist[2][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(kk);
                atom_.type = buf[i].type;
                atom_.x[0] = buf[i].r[0];
                atom_.x[1] = buf[i].r[1];
                atom_.x[2] = buf[i].r[2];
            }
        } else {
            for (int i = 0; i < n; i++) {
                kk = recvlist[3][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(kk);
                atom_.type = buf[i].type;
                atom_.x[0] = buf[i].r[0];
                atom_.x[1] = buf[i].r[1];
                atom_.x[2] = buf[i].r[2];
            }
        }
    } else {
        if (direction == 0) {
            for (int i = 0; i < n; i++) {
                kk = recvlist[4][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(kk);
                atom_.type = buf[i].type;
                atom_.x[0] = buf[i].r[0];
                atom_.x[1] = buf[i].r[1];
                atom_.x[2] = buf[i].r[2];
            }
        } else {
            for (int i = 0; i < n; i++) {
                kk = recvlist[5][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(kk);
                atom_.type = buf[i].type;
                atom_.x[0] = buf[i].r[0];
                atom_.x[1] = buf[i].r[1];
                atom_.x[2] = buf[i].r[2];
            }
        }
    }
}

void pack::pack_rho(int n, AtomList &atom_list, double *buf, std::vector<_type_atom_id> &recvlist) {
    int j, m = 0;
    for (int i = 0; i < n; i++) {
        j = recvlist[i];
        AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
        buf[m++] = atom_.rho;
    }
}

void pack::unpack_rho(int d, int direction, AtomList &atom_list,
                      double *buf, std::vector<std::vector<_type_atom_id> > &sendlist) {
    int j, m = 0;
    if (d == 0) {
        if (direction == 0) {
            for (int i = 0; i < sendlist[1].size(); i++) {
                j = sendlist[1][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[0].size(); i++) {
                j = sendlist[0][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        }
    } else if (d == 1) {
        if (direction == 0) {
            for (int i = 0; i < sendlist[3].size(); i++) {
                j = sendlist[3][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[2].size(); i++) {
                j = sendlist[2][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        }
    } else {
        if (direction == 0) {
            for (int i = 0; i < sendlist[5].size(); i++) {
                j = sendlist[5][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[4].size(); i++) {
                j = sendlist[4][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        }
    }
}

void pack::pack_df(AtomList &atom_list, double *buf,
//                   int &nlocalinter,
//                   int &nghostinter,
//                   std::vector<_type_atom_id> &idinter,
//                   std::vector<_type_atom_type> &typeinter,
//                   std::vector<std::vector<double>> &xinter,
//                   std::vector<std::vector<double>> &vinter,
//                   std::vector<std::vector<double>> &finter,
//                   std::vector<double> &rhointer,
                   std::vector<double> &dfinter,
                   std::vector<_type_atom_id> &sendlist, std::vector<int> &intersendlist) {
    int j, m = 0;
    int n = sendlist.size();
    for (int i = 0; i < n; i++) {
        j = sendlist[i];
        AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
        buf[m++] = atom_.df;
    }
    n = intersendlist.size();
    for (int i = 0; i < n; i++) {
        j = intersendlist[i];
        buf[m++] = dfinter[j];
    }
}

void pack::unpack_df(int n, AtomList &atom_list, double *buf,
//                     int &nlocalinter,
//                     int &nghostinter,
//                     std::vector<_type_atom_id> &idinter,
//                     std::vector<_type_atom_type> &typeinter,
//                     std::vector<std::vector<double>> &xinter,
//                     std::vector<std::vector<double>> &vinter,
//                     std::vector<std::vector<double>> &finter,
//                     std::vector<double> &rhointer,
                     std::vector<double> &dfinter,
                     std::vector<_type_atom_id> &recvlist, std::vector<int> &interrecvlist) {
    long kk;
    int m = 0;
    if (n != (recvlist.size() + interrecvlist.size())) {
        printf("wrong number of dfembed recv!!!");
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
    n = recvlist.size();
    for (int i = 0; i < n; i++) {
        kk = recvlist[i];
        AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(kk);
        atom_.df = buf[m++];
    }
    n = interrecvlist.size();
    for (int i = 0; i < n; i++) {
        kk = interrecvlist[i];
        if (kk == -1) {
            m++;
            continue;
        }
        dfinter[kk] = buf[m++];
    }
}


void pack::pack_force(int n, AtomList &atom_list, double *buf,
                      std::vector<_type_atom_id> &recvlist) {
    int j, m = 0;
    for (int i = 0; i < n; i++) {
        j = recvlist[i];
        AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
        buf[m++] = atom_.f[0];
        buf[m++] = atom_.f[1];
        buf[m++] = atom_.f[2];
    }
}

void pack::unpack_force(int d, int direction, AtomList &atom_list, double *buf,
                        std::vector<std::vector<_type_atom_id> > &sendlist) {
    int j, m = 0;
    if (d == 0) {
        if (direction == 0) {
            for (int i = 0; i < sendlist[1].size(); i++) {
                j = sendlist[1][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.f[0] += buf[m++];
                atom_.f[1] += buf[m++];
                atom_.f[2] += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[0].size(); i++) {
                j = sendlist[0][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.f[0] += buf[m++];
                atom_.f[1] += buf[m++];
                atom_.f[2] += buf[m++];
            }
        }
    } else if (d == 1) {
        if (direction == 0) {
            for (int i = 0; i < sendlist[3].size(); i++) {
                j = sendlist[3][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.f[0] += buf[m++];
                atom_.f[1] += buf[m++];
                atom_.f[2] += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[2].size(); i++) {
                j = sendlist[2][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.f[0] += buf[m++];
                atom_.f[1] += buf[m++];
                atom_.f[2] += buf[m++];
            }
        }
    } else {
        if (direction == 0) {
            for (int i = 0; i < sendlist[5].size(); i++) {
                j = sendlist[5][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.f[0] += buf[m++];
                atom_.f[1] += buf[m++];
                atom_.f[2] += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[4].size(); i++) {
                j = sendlist[4][i];
                AtomElement &atom_ = atom_list.getAtomEleByLinearIndex(j);
                atom_.f[0] += buf[m++];
                atom_.f[1] += buf[m++];
                atom_.f[2] += buf[m++];
            }
        }
    }
}