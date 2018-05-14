#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <logs/logs.h>
#include "atom.h"
#include "toml_config.h"
#include "hardware_accelerate.hpp" // use hardware(eg.GPU, MIC,Sunway slave cores.) to achieve calculate accelerating.

atom::atom(Domain *domain, double latticeconst,
           double cutoffRadiusFactor, int seed) :
        p_domain(domain), _latticeconst(latticeconst),
        _cutoffRadius(cutoffRadiusFactor * latticeconst), _seed(seed) {

    _cutlattice = static_cast<int>(ceil(cutoffRadiusFactor));

    nlocalinter = 0;
    nghostinter = 0;

    numberoflattice =
            p_domain->getGhostLatticeSize(0) * p_domain->getGhostLatticeSize(1) * p_domain->getGhostLatticeSize(2);
    // printf("number:%d, %d, %d, %d, %d\n", numberoflattice,
    // p_domain->getSubBoxLatticeSize(0), p_domain->getSubBoxLatticeSize(1), p_domain->getSubBoxLatticeSize(2),
    // p_domain->getSubBoxLatticeSize(0)*p_domain->getSubBoxLatticeSize(1)*p_domain->getSubBoxLatticeSize(2));
    atom_list = new AtomList(p_domain->getGhostLatticeSize(0),
                             p_domain->getGhostLatticeSize(1),
                             p_domain->getGhostLatticeSize(2),
                             p_domain->getSubBoxLatticeCoordLower(0) - p_domain->getGhostLatticeCoordLower(0),
                             p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1),
                             p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2));

    calculateNeighbourIndices();

    if (isAccelerateSupport()) {
        accelerateInit(p_domain->getSubBoxLatticeCoordLower(0),
                       p_domain->getSubBoxLatticeCoordLower(1),
                       p_domain->getSubBoxLatticeCoordLower(2),
                       p_domain->getSubBoxLatticeSize(0),
                       p_domain->getSubBoxLatticeSize(1),
                       p_domain->getSubBoxLatticeSize(2),
                       p_domain->getGhostLatticeCoordLower(0),
                       p_domain->getGhostLatticeCoordLower(1),
                       p_domain->getGhostLatticeCoordLower(2),
                       p_domain->getGhostLatticeSize(0),
                       p_domain->getGhostLatticeSize(1),
                       p_domain->getGhostLatticeSize(2));
    }
}

atom::~atom() {
    delete atom_list;
}

void atom::calculateNeighbourIndices() {
    double x, y, z;
    int mark = 0;
    double cut_times_lattice = _cutoffRadius / _latticeconst;
    vector<long int>::iterator neighbourOffsetsIter;
    for (int zIndex = -_cutlattice; zIndex <= _cutlattice; zIndex++) { // loop for (2*_cutlattice + 1) times.
        for (int yIndex = -_cutlattice; yIndex <= _cutlattice; yIndex++) {
            for (int xIndex = -_cutlattice * 2; xIndex <= _cutlattice * 2; xIndex++) {
                // 体心
                z = (double) zIndex + (((double) (xIndex % 2)) / 2); // zIndex plus 1/2 (odd) or 0(even).
                y = (double) yIndex + (((double) (xIndex % 2)) / 2);
                x = (double) xIndex / 2;
                long int offset;
                double r = sqrt(x * x + y * y + z * z);
                if (r < (cut_times_lattice + 0.4)) { // todo 0.4?
                    offset = IndexOf3DIndex(xIndex, yIndex, zIndex);
                    if (offset > 0) {
                        NeighbourOffsets.push_back(offset);
                    }
                }

                // 晶格点
                z = (double) zIndex - (((double) (xIndex % 2)) / 2);
                y = (double) yIndex - (((double) (xIndex % 2)) / 2);
                x = (double) xIndex / 2;
                r = sqrt(x * x + y * y + z * z);
                if (r < (cut_times_lattice + 0.4)) {
                    offset = IndexOf3DIndex(xIndex, yIndex, zIndex);
                    if (offset > 0) {
                        for (neighbourOffsetsIter = NeighbourOffsets.begin();
                             neighbourOffsetsIter != NeighbourOffsets.end(); neighbourOffsetsIter++) {
                            if (*neighbourOffsetsIter == offset) {
                                mark = 1;
                            }
                        }
                        if (mark != 1) {
                            NeighbourOffsets.push_back(offset);
                        }
                        mark = 0;
                    }
                }
            }
        }
    }
}

long atom::IndexOf3DIndex(long int xIndex, long int yIndex, long int zIndex) const {
    return (zIndex * p_domain->getGhostLatticeSize(1) + yIndex) * p_domain->getGhostLatticeSize(0) + xIndex;
}

void atom::addAtom(_type_atom_id id, double rx, double ry, double rz, double vx, double vy, double vz) {
    int i;
    if ((rx >= p_domain->getMeasuredSubBoxLowerBounding(0)) &&
        (rx < p_domain->getMeasuredSubBoxUpperBounding(0)) &&
        (ry >= p_domain->getMeasuredSubBoxLowerBounding(1)) &&
        (ry < p_domain->getMeasuredSubBoxUpperBounding(1)) &&
        (rz >= p_domain->getMeasuredSubBoxLowerBounding(2)) &&
        (rz < p_domain->getMeasuredSubBoxUpperBounding(2))) {
        int lattice[3];
        lattice[0] = rx * 2 / _latticeconst + 0.5;
        lattice[1] = ry * 2 / _latticeconst + 0.5;
        lattice[2] = rz * 2 / _latticeconst + 0.5;
        lattice[1] = lattice[1] / 2;
        lattice[2] = lattice[2] / 2;
        lattice[0] -= p_domain->getGhostLatticeCoordLower(0);
        lattice[1] -= p_domain->getGhostLatticeCoordLower(1);
        lattice[2] -= p_domain->getGhostLatticeCoordLower(2);
        i = ((p_domain->getGhostLatticeSize(1)) * lattice[2] + lattice[1]) *
            (p_domain->getGhostLatticeSize(0)) + lattice[0];
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(i);
        atom_.id = id;
        atom_.x[0] = rx;
        atom_.x[1] = ry;
        atom_.x[2] = rz;
        atom_.v[0] = vx;
        atom_.v[1] = vy;
        atom_.v[2] = vz;
    }
}

int atom::decide() {
    nghostinter = 0;
    int nflag = 0;
    long kk = 0;
    double dist;
    double xtemp, ytemp, ztemp;
    int xstart = p_domain->getSubBoxLatticeCoordLower(0) - p_domain->getGhostLatticeCoordLower(0);
    int ystart = p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1);
    int zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);

    //对本地晶格点原子进行判断，看是否运动为间隙原子
    for (int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
        for (int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
            for (int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                kk = IndexOf3DIndex(i, j, k);
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                if (atom_.x[0] != COORDINATE_ATOM_OUT_BOX) {
                    xtemp = (i + p_domain->getGhostLatticeCoordLower(0)) * 0.5 * _latticeconst;
                    ytemp = (j + p_domain->getGhostLatticeCoordLower(1) + (i % 2) * 0.5) * _latticeconst;
                    ztemp = (k + p_domain->getGhostLatticeCoordLower(2) + (i % 2) * 0.5) * _latticeconst;
                    dist = (atom_.x[0] - xtemp) * (atom_.x[0] - xtemp);
                    dist += (atom_.x[1] - ytemp) * (atom_.x[1] - ytemp);
                    dist += (atom_.x[2] - ztemp) * (atom_.x[2] - ztemp);
                    if (dist > (pow(0.2 * _latticeconst, 2.0))) { /**超过距离则判断为间隙原子*/
                        if (xinter.size() > nlocalinter) {
                            if (idinter.size() > nlocalinter) {
                                idinter[nlocalinter] = atom_.id;
                            } else {
                                idinter.push_back(atom_.id);
                            }
                            typeinter[nlocalinter] = atom_.type;
                            xinter[nlocalinter][0] = atom_.x[0];
                            xinter[nlocalinter][1] = atom_.x[1];
                            xinter[nlocalinter][2] = atom_.x[2];

                            if (vinter.size() > nlocalinter) {
                                vinter[nlocalinter][0] = atom_.v[0];
                                vinter[nlocalinter][1] = atom_.v[1];
                                vinter[nlocalinter][2] = atom_.v[2];
                            } else {
                                vinter.resize(nlocalinter + 1, vector<double>(3));
                                vinter[nlocalinter][0] = atom_.v[0];
                                vinter[nlocalinter][1] = atom_.v[1];
                                vinter[nlocalinter][2] = atom_.v[2];
                            }
                            nlocalinter++;
                            finter.resize(nlocalinter, vector<double>(3));
                            rhointer.resize(nlocalinter);
                            dfinter.resize(nlocalinter);
                        } else {
                            idinter.push_back(atom_.id);
                            typeinter.push_back(atom_.type);
                            xinter.resize(nlocalinter + 1, vector<double>(3));
                            xinter[nlocalinter][0] = atom_.x[0];
                            xinter[nlocalinter][1] = atom_.x[1];
                            xinter[nlocalinter][2] = atom_.x[2];
                            vinter.resize(nlocalinter + 1, vector<double>(3));
                            vinter[nlocalinter][0] = atom_.v[0];
                            vinter[nlocalinter][1] = atom_.v[1];
                            vinter[nlocalinter][2] = atom_.v[2];
                            nlocalinter++;
                            finter.resize(nlocalinter, vector<double>(3));
                            rhointer.resize(nlocalinter);
                            dfinter.resize(nlocalinter);
                        }

                        atom_.x[0] = COORDINATE_ATOM_OUT_BOX;
                        atom_.x[1] = COORDINATE_ATOM_OUT_BOX;
                        atom_.x[2] = COORDINATE_ATOM_OUT_BOX;
                        atom_.v[0] = 0;
                        atom_.v[1] = 0;
                        atom_.v[2] = 0;
                        nflag = 1;
                    }
                }
            }
        }
    }

    // periodic boundary
    for (int i = 0; i < nlocalinter; i++) {
        if (xinter[i][0] < p_domain->getMeasuredGlobalBoxCoordLower(0)) {
            xinter[i][0] += p_domain->getMeasuredGlobalLength(0);
        } else if (xinter[i][0] >= p_domain->getMeasuredGlobalBoxCoordUpper(0)) {
            xinter[i][0] -= p_domain->getMeasuredGlobalLength(0);
        }
        if (xinter[i][1] < p_domain->getMeasuredGlobalBoxCoordLower(1)) {
            xinter[i][1] += p_domain->getMeasuredGlobalLength(1);
        } else if (xinter[i][1] >= p_domain->getMeasuredGlobalBoxCoordUpper(1)) {
            xinter[i][1] -= p_domain->getMeasuredGlobalLength(1);
        }
        if (xinter[i][2] < p_domain->getMeasuredGlobalBoxCoordLower(1)) {
            xinter[i][2] += p_domain->getMeasuredGlobalLength(2);
        } else if (xinter[i][2] >= p_domain->getMeasuredGlobalBoxCoordUpper(1)) {
            xinter[i][2] -= p_domain->getMeasuredGlobalLength(2);
        }
    }

    //判断，如果跑出晶格点的?佑峙芑鼐Ц竦悖蚍呕鼐Ц竦闶榇娲⑵湫畔?
    for (int i = 0; i < nlocalinter; i++) {
        int j, k, l;
        xtemp = xinter[i][0];
        ytemp = xinter[i][1];
        ztemp = xinter[i][2];
        j = xtemp * 2 / _latticeconst + 0.5;
        k = ytemp * 2 / _latticeconst + 0.5;
        l = ztemp * 2 / _latticeconst + 0.5;
        k = k / 2;
        l = l / 2;
        j -= p_domain->getGhostLatticeCoordLower(0);
        k -= p_domain->getGhostLatticeCoordLower(1);
        l -= p_domain->getGhostLatticeCoordLower(2);
        //判断是否在所表示晶格范围内
        if (j <= (p_domain->getSubBoxLatticeSize(0) + 2 * (ceil(_cutoffRadius / _latticeconst) + 1))
            && k <= (p_domain->getSubBoxLatticeSize(1) + (ceil(_cutoffRadius / _latticeconst) + 1))
            && l <= (p_domain->getSubBoxLatticeSize(2) + (ceil(_cutoffRadius / _latticeconst) + 1))) {
            j = IndexOf3DIndex(j, k, l);
            AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
            if (atom_.x[0] == COORDINATE_ATOM_OUT_BOX) {
                atom_.id = idinter[i];
                atom_.type = typeinter[i];
                atom_.x[0] = xinter[i][0];
                atom_.x[1] = xinter[i][1];
                atom_.x[2] = xinter[i][2];
                atom_.v[0] = vinter[i][0];
                atom_.v[1] = vinter[i][1];
                atom_.v[2] = vinter[i][2];

                idinter[i] = idinter[nlocalinter - 1];
                typeinter[i] = typeinter[nlocalinter - 1];
                xinter[i][0] = xinter[nlocalinter - 1][0];
                xinter[i][1] = xinter[nlocalinter - 1][1];
                xinter[i][2] = xinter[nlocalinter - 1][2];
                vinter[i][0] = vinter[nlocalinter - 1][0];
                vinter[i][1] = vinter[nlocalinter - 1][1];
                vinter[i][2] = vinter[nlocalinter - 1][2];

                i--;
                nlocalinter--;
            }
        }
    }
    return nflag;
}

void atom::clearForce() {
    for (int i = 0; i < numberoflattice; i++) {
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(i);
        atom_.f[0] = 0;
        atom_.f[1] = 0;
        atom_.f[2] = 0;
        atom_.rho = 0;
    }
    for (int i = 0; i < finter.size(); i++) {
        finter[i][0] = 0;
        finter[i][1] = 0;
        finter[i][2] = 0;
        rhointer[i] = 0;
    }
}

void atom::computeEam(eam *pot, Domain *domain, double &comm) {
    double starttime, stoptime;
    double xtemp, ytemp, ztemp;
    double delx, dely, delz;
    vector<long int>::iterator neighbourOffsetsIter;
    InterpolationObject *rho_spline = pot->rho;
    InterpolationObject *f_spline = pot->f;
    InterpolationObject *phi_spline = pot->phi;
    int n;
    double dist2;
    double r;
    double rhoTmp, dRho, dEmbed, dfEmbed, phiTmp, dPhi;
    int nr, m;
    double p;
    double (*spline)[7];
    double fpair;
    double recip, phi, phip, psip, z2, z2p;
    long kk;
    int xstart = p_domain->getSubBoxLatticeCoordLower(0) - p_domain->getGhostLatticeCoordLower(0);
    int ystart = p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1);
    int zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);

    // 本地晶格点上的原子计算电子云密度
    if (isAccelerateSupport()) {
        accelerateEamRhoCalc(&(rho_spline->n), atom_list, &_cutoffRadius,
                             &(rho_spline->invDx), rho_spline->values); // fixme
    } else { // calculate rho use cpu only.
        for (int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
            for (int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
                for (int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                    kk = IndexOf3DIndex(i, j, k);
                    AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                    xtemp = atom_.x[0];
                    ytemp = atom_.x[1];
                    ztemp = atom_.x[2];
                    if (xtemp != COORDINATE_ATOM_OUT_BOX) {
                        //对晶格点邻居原子遍历
                        for (neighbourOffsetsIter = NeighbourOffsets.begin();
                             neighbourOffsetsIter != NeighbourOffsets.end(); neighbourOffsetsIter++) {
                            n = (kk + *neighbourOffsetsIter);
                            AtomElement &atom_n = atom_list->getAtomEleByLinearIndex(n);
                            delx = xtemp - atom_n.x[0];
                            dely = ytemp - atom_n.x[1];
                            delz = ztemp - atom_n.x[2];
                            dist2 = delx * delx + dely * dely + delz * delz;
                            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                                r = sqrt(dist2);
                                nr = rho_spline->n;
                                p = r * rho_spline->invDx + 1.0;
                                m = static_cast<int> (p);
                                m = std::max(1, std::min(m, (nr - 1)));
                                p -= m;
                                p = std::min(p, 1.0);
                                spline = rho_spline->spline;
                                rhoTmp = ((spline[m][3] * p + spline[m][4]) * p + spline[m][5]) * p +
                                         spline[m][6];

                                atom_.rho += rhoTmp;
                                atom_n.rho += rhoTmp;
                            }
                        }
                    }
                }
            }
        }
    }

    //间隙原子电子云密度
    int j, k, l;
    for (int i = 0; i < nlocalinter; i++) {
        xtemp = xinter[i][0];
        ytemp = xinter[i][1];
        ztemp = xinter[i][2];
        j = xtemp * 2 / _latticeconst + 0.5;
        k = ytemp * 2 / _latticeconst + 0.5;
        l = ztemp * 2 / _latticeconst + 0.5;
        k = k / 2;
        l = l / 2;
        j -= p_domain->getGhostLatticeCoordLower(0);
        k -= p_domain->getGhostLatticeCoordLower(1);
        l -= p_domain->getGhostLatticeCoordLower(2);
        j = IndexOf3DIndex(j, k, l);
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);

        delx = xtemp - atom_.x[0];
        dely = ytemp - atom_.x[1];
        delz = ztemp - atom_.x[2];
        dist2 = delx * delx + dely * dely + delz * delz;
        if (dist2 < (_cutoffRadius * _cutoffRadius)) {
            r = sqrt(dist2);
            nr = rho_spline->n;
            p = r * rho_spline->invDx + 1.0;
            m = static_cast<int> (p);
            m = std::max(1, std::min(m, (nr - 1)));
            p -= m;
            p = std::min(p, 1.0);
            spline = rho_spline->spline;
            rhoTmp = ((spline[m][3] * p + spline[m][4]) * p + spline[m][5]) * p + spline[m][6];

            rhointer[i] += rhoTmp;
            atom_.rho += rhoTmp;
        }
        for (neighbourOffsetsIter = NeighbourOffsets.begin();
             neighbourOffsetsIter != NeighbourOffsets.end(); neighbourOffsetsIter++) {
            //计算间隙原子的所有邻居
            n = (j + *neighbourOffsetsIter);
            AtomElement &atom_n = atom_list->getAtomEleByLinearIndex(n);
            delx = xtemp - atom_n.x[0];
            dely = ytemp - atom_n.x[1];
            delz = ztemp - atom_n.x[2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                r = sqrt(dist2);
                nr = rho_spline->n;
                p = r * rho_spline->invDx + 1.0;
                m = static_cast<int> (p);
                m = std::max(1, std::min(m, (nr - 1)));
                p -= m;
                p = std::min(p, 1.0);
                spline = rho_spline->spline;
                rhoTmp = ((spline[m][3] * p + spline[m][4]) * p + spline[m][5]) * p + spline[m][6];

                rhointer[i] += rhoTmp;
                atom_n.rho += rhoTmp;
            }
            n = (j - *neighbourOffsetsIter);
            AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(n);
            delx = xtemp - atom_n.x[0];
            dely = ytemp - atom_n.x[1];
            delz = ztemp - atom_n.x[2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                r = sqrt(dist2);
                nr = rho_spline->n;
                p = r * rho_spline->invDx + 1.0;
                m = static_cast<int> (p);
                m = std::max(1, std::min(m, (nr - 1)));
                p -= m;
                p = std::min(p, 1.0);
                spline = rho_spline->spline;
                rhoTmp = ((spline[m][3] * p + spline[m][4]) * p + spline[m][5]) * p + spline[m][6];

                rhointer[i] += rhoTmp;
                atom_n.rho += rhoTmp;
            }
        }
        //对间隙原子遍历
        for (int k = i + 1; k < (nghostinter + nlocalinter); k++) {
            delx = xtemp - xinter[k][0];
            dely = ytemp - xinter[k][1];
            delz = ztemp - xinter[k][2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                r = sqrt(dist2);
                nr = rho_spline->n;
                p = r * rho_spline->invDx + 1.0;
                m = static_cast<int> (p);
                m = std::max(1, std::min(m, (nr - 1)));
                p -= m;
                p = std::min(p, 1.0);
                spline = rho_spline->spline;
                rhoTmp = ((spline[m][3] * p + spline[m][4]) * p + spline[m][5]) * p + spline[m][6];

                rhointer[i] += rhoTmp;
                rhointer[k] += rhoTmp;
            }
        }
        //计算间隙原子嵌入能导数
        nr = f_spline->n;
        p = rhointer[i] * f_spline->invDx + 1.0;
        m = static_cast<int> (p);
        m = std::max(1, std::min(m, (nr - 1)));
        p -= m;
        p = std::min(p, 1.0);
        spline = f_spline->spline;
        dfEmbed = (spline[m][0] * p + spline[m][1]) * p + spline[m][2];
        dfinter[i] = dfEmbed;
    }

//    ofstream outfile;
    /* char tmp[20];
    sprintf(tmp, "rho.atom");
    outfile.open(tmp);
    for(int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++){
            for(int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++){
                    for(int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++){
                            kk = IndexOf3DIndex( i, j, k);
                             AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                            if(atom_.x[0] != COORDINATE_ATOM_OUT_BOX)
                                    outfile << atom_.rho << std::endl;
                    }
            }
    }
for(int i = 0; i < rho_spline->n; i++){ // 1.todo remove start.
    outfile << rho_spline->spline[i][6] << std::endl;
} // 1. todo remove end.
    outfile.close();*/

    //发送电子云密度
    starttime = MPI_Wtime();
    domain->sendrho(this);
    stoptime = MPI_Wtime();
    comm = stoptime - starttime;

    /*sprintf(tmp, "rho2.atom");
    outfile;
    outfile.open(tmp);
    for(int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++){
            for(int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++){
                    for(int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++){
                            kk = IndexOf3DIndex( i, j, k);
                             AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                            if(atom_.x[0] != COORDINATE_ATOM_OUT_BOX)
                                    outfile << atom_.rho << std::endl;
                    }
            }
    }
    outfile.close();*/

    //本地晶格点计算嵌入能导数
    if (isAccelerateSupport()) {
//        std::cout << "df\n";
        accelerateEamDfCalc(&(f_spline->n), atom_list, &_cutoffRadius,
                            &(f_spline->invDx), f_spline->values);
    } else {
        for (int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
            for (int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
                for (int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                    kk = IndexOf3DIndex(i, j, k);
                    AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                    nr = f_spline->n;
                    p = atom_.rho * f_spline->invDx + 1.0;
                    m = static_cast<int> (p);
                    m = std::max(1, std::min(m, (nr - 1)));
                    p -= m;
                    p = std::min(p, 1.0);
                    spline = f_spline->spline;
                    dfEmbed = (spline[m][0] * p + spline[m][1]) * p + spline[m][2];
                    atom_.df = dfEmbed;
                }
            }
        }
    }

    /*sprintf(tmp, "df.atom");
    outfile.open(tmp);
    for(int i = 0; i < f_spline->n; i++){
        outfile << i << " " << f_spline->spline[i][6] << std::endl;
    }
    outfile.close();*/

    //发送嵌入能导数
    starttime = MPI_Wtime();
    domain->sendDfEmbed(this);
    stoptime = MPI_Wtime();
    comm += stoptime - starttime;

    if (isAccelerateSupport()) {
//        std::cout << "f\n";
        accelerateEamForceCalc(&(phi_spline->n), atom_list, &_cutoffRadius,
                               &(phi_spline->invDx), phi_spline->values, rho_spline->values);
    } else {
        /*sprintf(tmp, "f.atom");
        outfile.open(tmp);
        for(int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++){
                for(int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++){
                        for(int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++){
                                kk = IndexOf3DIndex( i, j, k);
                                 AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                                if(atom_.x[0] != COORDINATE_ATOM_OUT_BOX)
                                        outfile << f[kk*3] << std::endl;
                        }
                }
        }
        outfile.close();*/

        for (int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
            for (int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
                for (int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                    kk = IndexOf3DIndex(i, j, k);
                    AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                    xtemp = atom_.x[0];
                    ytemp = atom_.x[1];
                    ztemp = atom_.x[2];
                    if (xtemp != COORDINATE_ATOM_OUT_BOX) {
                        //对晶格点邻居原子遍历
                        for (neighbourOffsetsIter = NeighbourOffsets.begin();
                             neighbourOffsetsIter != NeighbourOffsets.end(); neighbourOffsetsIter++) {
                            n = (kk + *neighbourOffsetsIter);
                            AtomElement &atom_n = atom_list->getAtomEleByLinearIndex(n);
                            delx = xtemp - atom_n.x[0];
                            dely = ytemp - atom_n.x[1];
                            delz = ztemp - atom_n.x[2];
                            dist2 = delx * delx + dely * dely + delz * delz;
                            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                                r = sqrt(dist2);
                                nr = phi_spline->n;
                                p = r * phi_spline->invDx + 1.0;
                                m = static_cast<int> (p);
                                m = std::max(1, std::min(m, (nr - 1)));
                                p -= m;
                                p = std::min(p, 1.0);
                                spline = phi_spline->spline;
                                phiTmp = ((spline[m][3] * p + spline[m][4]) * p + spline[m][5]) * p +
                                         spline[m][6];
                                dPhi = (spline[m][0] * p + spline[m][1]) * p + spline[m][2];
                                spline = rho_spline->spline;
                                dRho = (spline[m][0] * p + spline[m][1]) * p + spline[m][2];

                                z2 = phiTmp;
                                z2p = dPhi;
                                recip = 1.0 / r;
                                phi = z2 * recip;
                                phip = z2p * recip - phi * recip;
                                psip = (atom_.df + atom_n.df) * dRho + phip;
                                fpair = -psip * recip;

                                atom_.f[0] += delx * fpair;
                                atom_.f[1] += dely * fpair;
                                atom_.f[2] += delz * fpair;

                                atom_n.f[0] -= delx * fpair;
                                atom_n.f[1] -= dely * fpair;
                                atom_n.f[2] -= delz * fpair;
                            }
                        }
                    }
                }
            }
        }
    } // end of if-isAccelerateSupport.

    /*sprintf(tmp, "f.atom");  // 2.todo remove start.
      outfile.open(tmp);
      for(int i = 0; i < phi_spline->n; i++){
         outfile << i << " " << phi_spline->spline[i][6] << std::endl;
      }
      outfile.close();*/ // 2.todo remove end.

    //间隙原子计算嵌入能和对势带来的力
    for (int i = 0; i < nlocalinter; i++) {
        int j, k, l;
        xtemp = xinter[i][0];
        ytemp = xinter[i][1];
        ztemp = xinter[i][2];
        j = xtemp * 2 / _latticeconst + 0.5;
        k = ytemp * 2 / _latticeconst + 0.5;
        l = ztemp * 2 / _latticeconst + 0.5;
        k = k / 2;
        l = l / 2;
        j -= p_domain->getGhostLatticeCoordLower(0);
        k -= p_domain->getGhostLatticeCoordLower(1);
        l -= p_domain->getGhostLatticeCoordLower(2);
        j = IndexOf3DIndex(j, k, l);
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);

        delx = xtemp - atom_.x[0];
        dely = ytemp - atom_.x[1];
        delz = ztemp - atom_.x[2];
        dist2 = delx * delx + dely * dely + delz * delz;
        if (dist2 < (_cutoffRadius * _cutoffRadius)) {
            r = sqrt(dist2);
            nr = phi_spline->n;
            p = r * phi_spline->invDx + 1.0;
            m = static_cast<int> (p);
            m = std::max(1, std::min(m, (nr - 1)));
            p -= m;
            p = std::min(p, 1.0);
            spline = phi_spline->spline;
            phiTmp = ((spline[m][3] * p + spline[m][4]) * p + spline[m][5]) * p + spline[m][6];
            dPhi = (spline[m][0] * p + spline[m][1]) * p + spline[m][2];
            spline = rho_spline->spline;
            dRho = (spline[m][0] * p + spline[m][1]) * p + spline[m][2];

            z2 = phiTmp;
            z2p = dPhi;
            recip = 1.0 / r;
            phi = z2 * recip;
            phip = z2p * recip - phi * recip;
            psip = (dfinter[i] + atom_.df) * dRho + phip;
            fpair = -psip * recip;

            finter[i][0] += delx * fpair;
            finter[i][1] += dely * fpair;
            finter[i][2] += delz * fpair;

            atom_.f[0] -= delx * fpair;
            atom_.f[1] -= dely * fpair;
            atom_.f[2] -= delz * fpair;
        }
        for (neighbourOffsetsIter = NeighbourOffsets.begin();
             neighbourOffsetsIter != NeighbourOffsets.end(); neighbourOffsetsIter++) {
            //计算间隙原子的所有邻居
            n = (j + *neighbourOffsetsIter);
            AtomElement &atom_n = atom_list->getAtomEleByLinearIndex(n);
            delx = xtemp - atom_n.x[0];
            dely = ytemp - atom_n.x[1];
            delz = ztemp - atom_n.x[2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                r = sqrt(dist2);
                nr = phi_spline->n;
                p = r * phi_spline->invDx + 1.0;
                m = static_cast<int> (p);
                m = std::max(1, std::min(m, (nr - 1)));
                p -= m;
                p = std::min(p, 1.0);
                spline = phi_spline->spline;
                phiTmp = ((spline[m][3] * p + spline[m][4]) * p + spline[m][5]) * p + spline[m][6];
                dPhi = (spline[m][0] * p + spline[m][1]) * p + spline[m][2];
                spline = rho_spline->spline;
                dRho = (spline[m][0] * p + spline[m][1]) * p + spline[m][2];

                z2 = phiTmp;
                z2p = dPhi;
                recip = 1.0 / r;
                phi = z2 * recip;
                phip = z2p * recip - phi * recip;
                psip = (dfinter[i] + atom_n.df) * dRho + phip;
                fpair = -psip * recip;

                finter[i][0] += delx * fpair;
                finter[i][1] += dely * fpair;
                finter[i][2] += delz * fpair;

                atom_n.f[0] -= delx * fpair;
                atom_n.f[1] -= dely * fpair;
                atom_n.f[2] -= delz * fpair;
            }
            n = (j - *neighbourOffsetsIter);
            AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(n);

            delx = xtemp - atom_.x[0];
            dely = ytemp - atom_.x[1];
            delz = ztemp - atom_.x[2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                r = sqrt(dist2);
                nr = phi_spline->n;
                p = r * phi_spline->invDx + 1.0;
                m = static_cast<int> (p);
                m = std::max(1, std::min(m, (nr - 1)));
                p -= m;
                p = std::min(p, 1.0);
                spline = phi_spline->spline;
                phiTmp = ((spline[m][3] * p + spline[m][4]) * p + spline[m][5]) * p + spline[m][6];
                dPhi = (spline[m][0] * p + spline[m][1]) * p + spline[m][2];
                spline = rho_spline->spline;
                dRho = (spline[m][0] * p + spline[m][1]) * p + spline[m][2];

                z2 = phiTmp;
                z2p = dPhi;
                recip = 1.0 / r;
                phi = z2 * recip;
                phip = z2p * recip - phi * recip;
                psip = (dfinter[i] + atom_.df) * dRho + phip;
                fpair = -psip * recip;

                finter[i][0] += delx * fpair;
                finter[i][1] += dely * fpair;
                finter[i][2] += delz * fpair;

                atom_.f[0] -= delx * fpair;
                atom_.f[1] -= dely * fpair;
                atom_.f[2] -= delz * fpair;
            }
        }
        //对间隙原子遍历
        for (int k = i + 1; k < (nghostinter + nlocalinter); k++) {
            delx = xtemp - xinter[k][0];
            dely = ytemp - xinter[k][1];
            delz = ztemp - xinter[k][2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                r = sqrt(dist2);
                nr = phi_spline->n;
                p = r * phi_spline->invDx + 1.0;
                m = static_cast<int> (p);
                m = std::max(1, std::min(m, (nr - 1)));
                p -= m;
                p = std::min(p, 1.0);
                spline = phi_spline->spline;
                phiTmp = ((spline[m][3] * p + spline[m][4]) * p + spline[m][5]) * p + spline[m][6];
                dPhi = (spline[m][0] * p + spline[m][1]) * p + spline[m][2];
                spline = rho_spline->spline;
                dRho = (spline[m][0] * p + spline[m][1]) * p + spline[m][2];

                z2 = phiTmp;
                z2p = dPhi;
                recip = 1.0 / r;
                phi = z2 * recip;
                phip = z2p * recip - phi * recip;
                psip = (dfinter[i] + dfinter[k]) * dRho + phip;
                fpair = -psip * recip;

                finter[i][0] += delx * fpair;
                finter[i][1] += dely * fpair;
                finter[i][2] += delz * fpair;

                finter[k][0] -= delx * fpair;
                finter[k][1] -= dely * fpair;
                finter[k][2] -= delz * fpair;
            }
        }
    }
}

unsigned long atom::getinteridsendsize() {
    return interbuf.size();
}

void atom::getatomx(int direction, vector<vector<int> > &sendlist) {
    int i;
    if (direction == 0) {
        //找到要发送到邻居进程的区域
        int xstart = p_domain->getSubBoxLatticeCoordLower(0) - p_domain->getGhostLatticeCoordLower(0);
        int ystart = p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1);
        int zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);
        int xstop = xstart + (_cutlattice) * 2;
        int ystop = ystart + p_domain->getSubBoxLatticeSize(1);
        int zstop = zstart + p_domain->getSubBoxLatticeSize(2);

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[0].push_back(i);
                }
            }
        }
    } else {
        //找到要发送到邻居进程的区域
        int xstart = p_domain->getSubBoxLatticeCoordLower(0) - p_domain->getGhostLatticeCoordLower(0) +
                     p_domain->getSubBoxLatticeSize(0) - ((_cutlattice) * 2);
        int ystart = p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1);
        int zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);
        int xstop = p_domain->getSubBoxLatticeCoordLower(0) - p_domain->getGhostLatticeCoordLower(0) +
                    p_domain->getSubBoxLatticeSize(0);
        int ystop = ystart + p_domain->getSubBoxLatticeSize(1);
        int zstop = zstart + p_domain->getSubBoxLatticeSize(2);

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[1].push_back(i);
                }
            }
        }
    }
}

void atom::getatomy(int direction, vector<vector<int> > &sendlist) {
    int i;
    if (direction == 0) {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1);
        int zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);
        int xstop = p_domain->getGhostLatticeSize(0);
        int ystop = ystart + _cutlattice;
        int zstop = zstart + p_domain->getSubBoxLatticeSize(2);

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[2].push_back(i);
                }
            }
        }
    } else {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart =
                p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1) +
                p_domain->getSubBoxLatticeSize(1) - (_cutlattice);
        int zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);
        int xstop = p_domain->getGhostLatticeSize(0);
        int ystop = p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1) +
                    p_domain->getSubBoxLatticeSize(1);
        int zstop = zstart + p_domain->getSubBoxLatticeSize(2);

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[3].push_back(i);
                }
            }
        }
    }
}

void atom::getatomz(int direction, vector<vector<int> > &sendlist) {
    int i;
    if (direction == 0) {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = 0;
        int zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);
        int xstop = p_domain->getGhostLatticeSize(0);
        int ystop = p_domain->getGhostLatticeSize(1);
        int zstop = zstart + _cutlattice;

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[4].push_back(i);
                }
            }
        }
    } else {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = 0;
        int zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2) +
                     p_domain->getSubBoxLatticeSize(2) - (_cutlattice);
        int xstop = p_domain->getGhostLatticeSize(0);
        int ystop = p_domain->getGhostLatticeSize(1);
        int zstop = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2) +
                    p_domain->getSubBoxLatticeSize(2);

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = IndexOf3DIndex(ix, iy, iz);
                    sendlist[5].push_back(i);
                }
            }
        }
    }
}

void atom::getIntertosend(int d, int direction, double ghostlengh, vector<int> &sendlist) {
    double low, high;
    if (d == 0) {
        if (direction == 0) {
            low = p_domain->getMeasuredSubBoxLowerBounding(0);
            high = p_domain->getMeasuredSubBoxLowerBounding(0) + ghostlengh;
            for (int i = 0; i < nlocalinter; i++) {
                if (xinter[i][0] < high && xinter[i][0] >= low) {
                    sendlist.push_back(i);
                }
            }
        } else {
            low = p_domain->getMeasuredSubBoxUpperBounding(0) - ghostlengh;
            high = p_domain->getMeasuredSubBoxUpperBounding(0);
            for (int i = 0; i < nlocalinter; i++) {
                if (xinter[i][0] <= high && xinter[i][0] > low) {
                    sendlist.push_back(i);
                }
            }
        }
    } else if (d == 1) {
        if (direction == 0) {
            low = p_domain->getMeasuredSubBoxLowerBounding(1);
            high = p_domain->getMeasuredSubBoxLowerBounding(1) + ghostlengh;
            for (int i = 0; i < nlocalinter + nghostinter; i++) {
                if (xinter[i][1] < high && xinter[i][1] >= low) {
                    sendlist.push_back(i);
                }
            }
        } else {
            low = p_domain->getMeasuredSubBoxUpperBounding(1) - ghostlengh;
            high = p_domain->getMeasuredSubBoxUpperBounding(1);
            for (int i = 0; i < nlocalinter + nghostinter; i++) {
                if (xinter[i][1] <= high && xinter[i][1] > low) {
                    sendlist.push_back(i);
                }
            }
        }
    } else {
        if (direction == 0) {
            low = p_domain->getMeasuredSubBoxLowerBounding(2);
            high = p_domain->getMeasuredSubBoxLowerBounding(2) + ghostlengh;
            for (int i = 0; i < nlocalinter + nghostinter; i++) {
                if (xinter[i][2] < high && xinter[i][2] >= low) {
                    sendlist.push_back(i);
                }
            }
        } else {
            low = p_domain->getMeasuredSubBoxUpperBounding(2) - ghostlengh;
            high = p_domain->getMeasuredSubBoxUpperBounding(2);
            for (int i = 0; i < nlocalinter + nghostinter; i++) {
                if (xinter[i][2] <= high && xinter[i][2] > low) {
                    sendlist.push_back(i);
                }
            }
        }
    }
}

int atom::getintersendnum(int dimension, int direction) {
    interbuf.clear();
    for (int i = 0; i < nlocalinter; i++) {
        if (direction == 0) {
            if (xinter[i][dimension] < p_domain->getMeasuredSubBoxLowerBounding(dimension)) {
                interbuf.push_back(i);
            }
        } else {
            if (xinter[i][dimension] >= p_domain->getMeasuredSubBoxUpperBounding(dimension)) {
                interbuf.push_back(i);
            }
        }
    }
    return interbuf.size();
}

void atom::pack_intersend(particledata *buf) {
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

void atom::unpack_interrecv(int d, int n, particledata *buf) {
    vector<double> xtemp(3);
    vector<double> vtemp(3);
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
        if (xtemp[d] >= p_domain->getMeasuredSubBoxLowerBounding(d) &&
            xtemp[d] < p_domain->getMeasuredSubBoxUpperBounding(d)) {
            if (nlocalinter == xinter.size()) {
                idinter.push_back(id);
                typeinter.push_back(type);
                xinter.push_back(xtemp);
                vinter.push_back(vtemp);
                nlocalinter++;
                finter.resize(nlocalinter, vector<double>(3));
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
                finter.resize(nlocalinter, vector<double>(3));
                rhointer.resize(nlocalinter);
                dfinter.resize(nlocalinter);
            }
        }
    }
}

void atom::pack_bordersend(int dimension, int n, vector<int> &sendlist, LatParticleData *buf, double shift) {
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

void atom::unpack_borderrecv(int n, LatParticleData *buf, vector<int> &recvlist) {
    int type;
    vector<double> xtemp(3);
    for (int i = 0; i < n; i++) {
        type = buf[i].type;
        xtemp[0] = buf[i].r[0];
        xtemp[1] = buf[i].r[1];
        xtemp[2] = buf[i].r[2];
        if (xtemp[0] >= p_domain->getMeasuredGhostLowerBounding(0) &&
            xtemp[0] < p_domain->getMeasuredGhostUpperBounding(0) &&
            xtemp[1] >= p_domain->getMeasuredGhostLowerBounding(1) &&
            xtemp[1] < p_domain->getMeasuredGhostUpperBounding(1) &&
            xtemp[2] >= p_domain->getMeasuredGhostLowerBounding(2) &&
            xtemp[2] < p_domain->getMeasuredGhostUpperBounding(2)) {
            if (xinter.size() == nlocalinter + nghostinter) {
                typeinter.push_back(type);
                xinter.push_back(xtemp);
                nghostinter++;
                recvlist[i] = nlocalinter + nghostinter - 1;
                finter.resize(nlocalinter + nghostinter, vector<double>(3));
                rhointer.resize(nlocalinter + nghostinter);
                dfinter.resize(nlocalinter + nghostinter);
            } else {
                typeinter[nlocalinter + nghostinter] = type;
                xinter[nlocalinter + nghostinter][0] = xtemp[0];
                xinter[nlocalinter + nghostinter][1] = xtemp[1];
                xinter[nlocalinter + nghostinter][2] = xtemp[2];
                nghostinter++;
                recvlist[i] = nlocalinter + nghostinter - 1;
                finter.resize(nlocalinter + nghostinter, vector<double>(3));
                rhointer.resize(nlocalinter + nghostinter);
                dfinter.resize(nlocalinter + nghostinter);
            }
        } else
            recvlist[i] = -1;
    }
}

void atom::pack_send(int dimension, int n, vector<int> &sendlist, LatParticleData *buf, double shift) {
    int j;
    if (dimension == 0) {
        for (int i = 0; i < n; i++) {
            j = sendlist[i];
            AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
            buf[i].type = atom_.type;
            buf[i].r[0] = atom_.x[0] + shift;
            buf[i].r[1] = atom_.x[1];
            buf[i].r[2] = atom_.x[2];
        }
    } else if (dimension == 1) {
        for (int i = 0; i < n; i++) {
            j = sendlist[i];
            AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
            buf[i].type = atom_.type;
            buf[i].r[0] = atom_.x[0];
            buf[i].r[1] = atom_.x[1] + shift;
            buf[i].r[2] = atom_.x[2];
        }
    } else {
        for (int i = 0; i < n; i++) {
            j = sendlist[i];
            AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
            buf[i].type = atom_.type;
            buf[i].r[0] = atom_.x[0];
            buf[i].r[1] = atom_.x[1];
            buf[i].r[2] = atom_.x[2] + shift;
        }
    }
}

void atom::unpack_recvfirst(int d, int direction, int n, LatParticleData *buf, vector<vector<int> > &recvlist) {
    int xstart, ystart, zstart;
    int xstop, ystop, zstop;
    long kk;
    int m = 0;
    if (d == 0) {
        if (direction == 0) {
            xstart = p_domain->getSubBoxLatticeCoordLower(0) - p_domain->getGhostLatticeCoordLower(0) +
                     p_domain->getSubBoxLatticeSize(0);
            xstop = p_domain->getGhostLatticeSize(0);
            ystart = p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1);
            ystop = ystart + p_domain->getSubBoxLatticeSize(1);
            zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);
            zstop = zstart + p_domain->getSubBoxLatticeSize(2);
            for (int k = zstart; k < zstop; k++) {
                for (int j = ystart; j < ystop; j++) {
                    for (int i = xstart; i < xstop; i++) {
                        kk = IndexOf3DIndex(i, j, k);
                        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                        atom_.type = buf[m].type;
                        atom_.x[0] = buf[m].r[0];
                        atom_.x[1] = buf[m].r[1];
                        atom_.x[2] = buf[m++].r[2];
                        recvlist[0].push_back(kk);
                    }
                }
            }
            if (n != recvlist[0].size()) // todo error.
                printf("wrong!!!\n");
        } else {
            xstart = 0;
            xstop = p_domain->getSubBoxLatticeCoordLower(0) - p_domain->getGhostLatticeCoordLower(0);
            ystart = p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1);
            ystop = ystart + p_domain->getSubBoxLatticeSize(1);
            zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);
            zstop = zstart + p_domain->getSubBoxLatticeSize(2);
            for (int k = zstart; k < zstop; k++) {
                for (int j = ystart; j < ystop; j++) {
                    for (int i = xstart; i < xstop; i++) {
                        kk = IndexOf3DIndex(i, j, k);
                        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                        atom_.type = buf[m].type;
                        atom_.x[0] = buf[m].r[0];
                        atom_.x[1] = buf[m].r[1];
                        atom_.x[2] = buf[m++].r[2];
                        recvlist[1].push_back(kk);
                    }
                }
            }
            if (n != recvlist[1].size()) // todo error handling in dataReuse feature.
                printf("wrong!!!\n");
        }
    } else if (d == 1) {
        if (direction == 0) {
            xstart = 0;
            xstop = p_domain->getGhostLatticeSize(0);
            ystart = p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1) +
                     p_domain->getSubBoxLatticeSize(1);
            ystop = p_domain->getGhostLatticeSize(1);
            zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);
            zstop = zstart + p_domain->getSubBoxLatticeSize(2);
            for (int k = zstart; k < zstop; k++) {
                for (int j = ystart; j < ystop; j++) {
                    for (int i = xstart; i < xstop; i++) {
                        kk = IndexOf3DIndex(i, j, k);
                        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                        atom_.type = buf[m].type;
                        atom_.x[0] = buf[m].r[0];
                        atom_.x[1] = buf[m].r[1];
                        atom_.x[2] = buf[m++].r[2];
                        recvlist[2].push_back(kk);
                    }
                }
            }
            if (n != recvlist[2].size()) // todo error handling in dataReuse feature.
                printf("wrong!!!\n");
        } else {
            xstart = 0;
            xstop = p_domain->getGhostLatticeSize(0);
            ystart = 0;
            ystop = p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1);
            zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);
            zstop = zstart + p_domain->getSubBoxLatticeSize(2);
            for (int k = zstart; k < zstop; k++) {
                for (int j = ystart; j < ystop; j++) {
                    for (int i = xstart; i < xstop; i++) {
                        kk = IndexOf3DIndex(i, j, k);
                        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                        atom_.type = buf[m].type;
                        atom_.x[0] = buf[m].r[0];
                        atom_.x[1] = buf[m].r[1];
                        atom_.x[2] = buf[m++].r[2];
                        recvlist[3].push_back(kk);
                    }
                }
            }
            if (n != recvlist[3].size()) // todo error handling in dataReuse feature.
                printf("wrong!!!\n");
        }
    } else {
        if (direction == 0) {
            xstart = 0;
            xstop = p_domain->getGhostLatticeSize(0);
            ystart = 0;
            ystop = p_domain->getGhostLatticeSize(1);
            zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2) +
                     p_domain->getSubBoxLatticeSize(2);
            zstop = p_domain->getGhostLatticeSize(2);
            for (int k = zstart; k < zstop; k++) {
                for (int j = ystart; j < ystop; j++) {
                    for (int i = xstart; i < xstop; i++) {
                        kk = IndexOf3DIndex(i, j, k);
                        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                        atom_.type = buf[m].type;
                        atom_.x[0] = buf[m].r[0];
                        atom_.x[1] = buf[m].r[1];
                        atom_.x[2] = buf[m++].r[2];
                        recvlist[4].push_back(kk);
                    }
                }
            }
            if (n != recvlist[4].size()) // todo error handling in dataReuse feature.
                printf("wrong!!!\n");
        } else {
            xstart = 0;
            xstop = p_domain->getGhostLatticeSize(0);
            ystart = 0;
            ystop = p_domain->getGhostLatticeSize(1);
            zstart = 0;
            zstop = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);
            for (int k = zstart; k < zstop; k++) {
                for (int j = ystart; j < ystop; j++) {
                    for (int i = xstart; i < xstop; i++) {
                        kk = IndexOf3DIndex(i, j, k);
                        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                        atom_.type = buf[m].type;
                        atom_.x[0] = buf[m].r[0];
                        atom_.x[1] = buf[m].r[1];
                        atom_.x[2] = buf[m++].r[2];
                        recvlist[5].push_back(kk);
                    }
                }
            }
            if (n != recvlist[5].size()) // todo error handling in dataReuse feature.
                printf("wrong!!!\n");
        }
    }
}

void atom::unpack_recv(int d, int direction, int n, LatParticleData *buf, vector<vector<int> > &recvlist) {
    long kk;
    if (d == 0) {
        if (direction == 0) {
            for (int i = 0; i < n; i++) {
                kk = recvlist[0][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                atom_.type = buf[i].type;
                atom_.x[0] = buf[i].r[0];
                atom_.x[1] = buf[i].r[1];
                atom_.x[2] = buf[i].r[2];
            }
        } else {
            for (int i = 0; i < n; i++) {
                kk = recvlist[1][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
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
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                atom_.type = buf[i].type;
                atom_.x[0] = buf[i].r[0];
                atom_.x[1] = buf[i].r[1];
                atom_.x[2] = buf[i].r[2];
            }
        } else {
            for (int i = 0; i < n; i++) {
                kk = recvlist[3][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
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
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                atom_.type = buf[i].type;
                atom_.x[0] = buf[i].r[0];
                atom_.x[1] = buf[i].r[1];
                atom_.x[2] = buf[i].r[2];
            }
        } else {
            for (int i = 0; i < n; i++) {
                kk = recvlist[5][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                atom_.type = buf[i].type;
                atom_.x[0] = buf[i].r[0];
                atom_.x[1] = buf[i].r[1];
                atom_.x[2] = buf[i].r[2];
            }
        }
    }
}

void atom::pack_rho(int n, vector<int> &recvlist, double *buf) {
    int j, m = 0;
    for (int i = 0; i < n; i++) {
        j = recvlist[i];
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
        buf[m++] = atom_.rho;
    }
}

void atom::unpack_rho(int d, int direction, double *buf, vector<vector<int> > &sendlist) {
    int j, m = 0;
    if (d == 0) {
        if (direction == 0) {
            for (int i = 0; i < sendlist[1].size(); i++) {
                j = sendlist[1][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[0].size(); i++) {
                j = sendlist[0][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        }
    } else if (d == 1) {
        if (direction == 0) {
            for (int i = 0; i < sendlist[3].size(); i++) {
                j = sendlist[3][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[2].size(); i++) {
                j = sendlist[2][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        }
    } else {
        if (direction == 0) {
            for (int i = 0; i < sendlist[5].size(); i++) {
                j = sendlist[5][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[4].size(); i++) {
                j = sendlist[4][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
                atom_.rho += buf[m++];
            }
        }
    }
}

void atom::pack_df(vector<int> &sendlist, vector<int> &intersendlist, double *buf) {
    int j, m = 0;
    int n = sendlist.size();
    for (int i = 0; i < n; i++) {
        j = sendlist[i];
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
        buf[m++] = atom_.df;
    }
    n = intersendlist.size();
    for (int i = 0; i < n; i++) {
        j = intersendlist[i];
        buf[m++] = dfinter[j];
    }
}

void atom::unpack_df(int n, double *buf, vector<int> &recvlist, vector<int> &interrecvlist) {
    long kk;
    int m = 0;
    if (n != (recvlist.size() + interrecvlist.size())) {
        printf("wrong number of dfembed recv!!!");
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
    n = recvlist.size();
    for (int i = 0; i < n; i++) {
        kk = recvlist[i];
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
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

void atom::pack_force(int n, vector<int> &recvlist, double *buf) {
    int j, m = 0;
    for (int i = 0; i < n; i++) {
        j = recvlist[i];
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
        buf[m++] = atom_.f[0];
        buf[m++] = atom_.f[1];
        buf[m++] = atom_.f[2];
    }
}

void atom::unpack_force(int d, int direction, double *buf, vector<vector<int> > &sendlist) {
    int j, m = 0;
    if (d == 0) {
        if (direction == 0) {
            for (int i = 0; i < sendlist[1].size(); i++) {
                j = sendlist[1][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
                atom_.f[0] += buf[m++];
                atom_.f[1] += buf[m++];
                atom_.f[2] += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[0].size(); i++) {
                j = sendlist[0][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
                atom_.f[0] += buf[m++];
                atom_.f[1] += buf[m++];
                atom_.f[2] += buf[m++];
            }
        }
    } else if (d == 1) {
        if (direction == 0) {
            for (int i = 0; i < sendlist[3].size(); i++) {
                j = sendlist[3][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
                atom_.f[0] += buf[m++];
                atom_.f[1] += buf[m++];
                atom_.f[2] += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[2].size(); i++) {
                j = sendlist[2][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
                atom_.f[0] += buf[m++];
                atom_.f[1] += buf[m++];
                atom_.f[2] += buf[m++];
            }
        }
    } else {
        if (direction == 0) {
            for (int i = 0; i < sendlist[5].size(); i++) {
                j = sendlist[5][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
                atom_.f[0] += buf[m++];
                atom_.f[1] += buf[m++];
                atom_.f[2] += buf[m++];
            }
        } else {
            for (int i = 0; i < sendlist[4].size(); i++) {
                j = sendlist[4][i];
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
                atom_.f[0] += buf[m++];
                atom_.f[1] += buf[m++];
                atom_.f[2] += buf[m++];
            }
        }
    }
}

void atom::computefirst(double dtInv2m, double dt) {
    long kk;
    int xstart = p_domain->getSubBoxLatticeCoordLower(0) - p_domain->getGhostLatticeCoordLower(0);
    int ystart = p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1);
    int zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);
    //本地晶格点上的原子求解运动方程第一步
    for (int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
        for (int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
            for (int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                kk = IndexOf3DIndex(i, j, k);
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                if (atom_.x[0] != COORDINATE_ATOM_OUT_BOX) {
                    atom_.v[0] = atom_.v[0] + dtInv2m * atom_.f[0];
                    atom_.x[0] += dt * atom_.v[0];
                    atom_.v[1] = atom_.v[1] + dtInv2m * atom_.f[1];
                    atom_.x[1] += dt * atom_.v[1];
                    atom_.v[2] = atom_.v[2] + dtInv2m * atom_.f[2];
                    atom_.x[2] += dt * atom_.v[2];
                }
            }
        }
    }
    //本地间隙原子求解运动方程第一步
    for (int i = 0; i < nlocalinter; i++) {
        for (unsigned short d = 0; d < 3; ++d) {
            vinter[i][d] = vinter[i][d] + dtInv2m * finter[i][d];
            xinter[i][d] += dt * vinter[i][d];
        }
    }
}

void atom::computesecond(double dtInv2m) {
    long kk;
    int xstart = p_domain->getSubBoxLatticeCoordLower(0) - p_domain->getGhostLatticeCoordLower(0);
    int ystart = p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1);
    int zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);
    //本地晶格点上的原子求解运动方程第二步
    for (int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
        for (int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
            for (int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                kk = IndexOf3DIndex(i, j, k);
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                for (unsigned short d = 0; d < 3; ++d) {
                    atom_.v[d] += dtInv2m * atom_.f[d]; // fixme hot!!
                }
            }
        }
    }
    //本地间隙原子求解运动方程第二步
    for (int i = 0; i < nlocalinter; i++) {
        for (unsigned short d = 0; d < 3; ++d) {
            vinter[i][d] = vinter[i][d] + dtInv2m * finter[i][d];
        }
    }
}

void atom::print_force() {
    char tmp[20];
    sprintf(tmp, "force.txt");
    ofstream outfile;
    outfile.open(tmp);
    long kk;
    int xstart = p_domain->getSubBoxLatticeCoordLower(0) - p_domain->getGhostLatticeCoordLower(0);
    int ystart = p_domain->getSubBoxLatticeCoordLower(1) - p_domain->getGhostLatticeCoordLower(1);
    int zstart = p_domain->getSubBoxLatticeCoordLower(2) - p_domain->getGhostLatticeCoordLower(2);
    std::cout << "print_force" << std::endl;
    for (int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
        for (int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
            for (int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                kk = IndexOf3DIndex(i, j, k);
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                outfile << atom_.f[0] << " " << atom_.f[1] << " " << atom_.f[2] << std::endl;
            }
        }
    }
    outfile.close();
}

void atom::setv(int lat[4], double collision_v[3]) {
    long kk;
    if ((lat[0] * 2) >= p_domain->getSubBoxLatticeCoordLower(0) &&
        (lat[0] * 2) < (p_domain->getSubBoxLatticeCoordLower(0) + p_domain->getSubBoxLatticeSize(0))
        && lat[1] >= p_domain->getSubBoxLatticeCoordLower(1) &&
        lat[1] < (p_domain->getSubBoxLatticeCoordLower(1) + p_domain->getSubBoxLatticeSize(1))
        && lat[2] >= p_domain->getSubBoxLatticeCoordLower(2) &&
        lat[2] < (p_domain->getSubBoxLatticeCoordLower(2) + p_domain->getSubBoxLatticeSize(2))) {
        kk = (IndexOf3DIndex(lat[0] * 2 - p_domain->getGhostLatticeCoordLower(0),
                             lat[1] - p_domain->getGhostLatticeCoordLower(1),
                             lat[2] - p_domain->getGhostLatticeCoordLower(2)) + lat[4]);
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
        atom_.v[0] += collision_v[0];
        atom_.v[1] += collision_v[1];
        atom_.v[2] += collision_v[2];
    }
}

int atom::getnlocalatom() {
    return (p_domain->getSubBoxLatticeSize(0) * p_domain->getSubBoxLatticeSize(1) *
            p_domain->getSubBoxLatticeSize(2));
}
