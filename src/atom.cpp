#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <logs/logs.h>
#include "atom.h"
#include "toml_config.h"
#include "hardware_accelerate.hpp" // use hardware(eg.GPU, MIC,Sunway slave cores.) to achieve calculate accelerating.
#include "potential/eam.h"

atom::atom(Domain *domain, double latticeconst,
           double cutoffRadiusFactor, int seed) :
        p_domain(domain), _latticeconst(latticeconst),
        _cutoffRadius(cutoffRadiusFactor * latticeconst), _seed(seed) {

    _cutlattice = static_cast<int>(ceil(cutoffRadiusFactor));

    numberoflattice = p_domain->getGhostExtLatticeSize(0) * p_domain->getGhostExtLatticeSize(1) *
                      p_domain->getGhostExtLatticeSize(2);
    // printf("number:%d, %d, %d, %d, %d\n", numberoflattice,
    // p_domain->getSubBoxLatticeSize(0), p_domain->getSubBoxLatticeSize(1), p_domain->getSubBoxLatticeSize(2),
    // p_domain->getSubBoxLatticeSize(0)*p_domain->getSubBoxLatticeSize(1)*p_domain->getSubBoxLatticeSize(2));
    atom_list = new AtomList(p_domain->getGhostExtLatticeSize(0),
                             p_domain->getGhostExtLatticeSize(1),
                             p_domain->getGhostExtLatticeSize(2),
                             p_domain->getSubBoxLatticeSize(0),
                             p_domain->getSubBoxLatticeSize(1),
                             p_domain->getSubBoxLatticeSize(2),
                             p_domain->getGhostLatticeSize(0),
                             p_domain->getGhostLatticeSize(1),
                             p_domain->getGhostLatticeSize(2));

    inter_atom_list = new InterAtomList();

    if (isAccelerateSupport()) {
        accelerateInit(p_domain->getGlobalSubBoxLatticeCoordLower(0),
                       p_domain->getGlobalSubBoxLatticeCoordLower(1),
                       p_domain->getGlobalSubBoxLatticeCoordLower(2),
                       p_domain->getSubBoxLatticeSize(0),
                       p_domain->getSubBoxLatticeSize(1),
                       p_domain->getSubBoxLatticeSize(2),
                       p_domain->getGlobalGhostLatticeCoordLower(0),
                       p_domain->getGlobalGhostLatticeCoordLower(1),
                       p_domain->getGlobalGhostLatticeCoordLower(2),
                       p_domain->getGhostExtLatticeSize(0),
                       p_domain->getGhostExtLatticeSize(1),
                       p_domain->getGhostExtLatticeSize(2));
    }
}

atom::~atom() {
    delete atom_list;
    delete inter_atom_list;
}

void atom::calculateNeighbourIndices() {
    double x, y, z;
    int mark = 0;
    double cut_times_lattice = _cutoffRadius / _latticeconst; // todo use cutoffRadiusFactor.
    std::vector<long int>::iterator neighbourOffsetsIter;
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
                    offset = atom_list->IndexOf3DIndex(xIndex, yIndex, zIndex);
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
                    offset = atom_list->IndexOf3DIndex(xIndex, yIndex, zIndex);
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
        lattice[0] -= p_domain->getGlobalGhostLatticeCoordLower(0);
        lattice[1] -= p_domain->getGlobalGhostLatticeCoordLower(1);
        lattice[2] -= p_domain->getGlobalGhostLatticeCoordLower(2);
        i = ((p_domain->getGhostExtLatticeSize(1)) * lattice[2] + lattice[1]) *
            (p_domain->getGhostExtLatticeSize(0)) + lattice[0];
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
    inter_atom_list->nghostinter = 0;
    int nflag = 0;
    long kk = 0;
    double dist;
    double xtemp, ytemp, ztemp;
//    int xstart = p_domain->getGhostLatticeSize(0);
//    int ystart = p_domain->getGhostLatticeSize(1);
//    int zstart = p_domain->getGhostLatticeSize(2);

    //对本地晶格点原子进行判断，看是否运动为间隙原子
    for (int k = 0; k < p_domain->getSubBoxLatticeSize(2); k++) {
        for (int j = 0; j < p_domain->getSubBoxLatticeSize(1); j++) {
            for (int i = 0; i < p_domain->getSubBoxLatticeSize(0); i++) {
//                kk = atom_list->IndexOf3DIndex(i, j, k);
                AtomElement &atom_ = atom_list->getAtomEleBySubBoxIndex(i, j, k); // todo long type
                if (!atom_.isInterElement()) {
                    xtemp = (i + p_domain->getGlobalSubBoxLatticeCoordLower(0)) * 0.5 * _latticeconst;
                    ytemp = (j + p_domain->getGlobalSubBoxLatticeCoordLower(1) + (i % 2) * 0.5) * _latticeconst;
                    ztemp = (k + p_domain->getGlobalSubBoxLatticeCoordLower(2) + (i % 2) * 0.5) * _latticeconst;
                    dist = (atom_.x[0] - xtemp) * (atom_.x[0] - xtemp);
                    dist += (atom_.x[1] - ytemp) * (atom_.x[1] - ytemp);
                    dist += (atom_.x[2] - ztemp) * (atom_.x[2] - ztemp);
                    if (dist > (pow(0.2 * _latticeconst, 2.0))) { /**超过距离则判断为间隙原子*/
                        if (inter_atom_list->xinter.size() > inter_atom_list->nlocalinter) {
                            if (inter_atom_list->idinter.size() > inter_atom_list->nlocalinter) {
                                inter_atom_list->idinter[inter_atom_list->nlocalinter] = atom_.id;
                            } else {
                                inter_atom_list->idinter.push_back(atom_.id);
                            }
                            inter_atom_list->typeinter[inter_atom_list->nlocalinter] = atom_.type; // todo type.
                            inter_atom_list->xinter[inter_atom_list->nlocalinter][0] = atom_.x[0];
                            inter_atom_list->xinter[inter_atom_list->nlocalinter][1] = atom_.x[1];
                            inter_atom_list->xinter[inter_atom_list->nlocalinter][2] = atom_.x[2];

                            if (inter_atom_list->vinter.size() > inter_atom_list->nlocalinter) {
                                inter_atom_list->vinter[inter_atom_list->nlocalinter][0] = atom_.v[0];
                                inter_atom_list->vinter[inter_atom_list->nlocalinter][1] = atom_.v[1];
                                inter_atom_list->vinter[inter_atom_list->nlocalinter][2] = atom_.v[2];
                            } else {
                                inter_atom_list->vinter.resize(inter_atom_list->nlocalinter + 1,
                                                               std::vector<double>(3));
                                inter_atom_list->vinter[inter_atom_list->nlocalinter][0] = atom_.v[0];
                                inter_atom_list->vinter[inter_atom_list->nlocalinter][1] = atom_.v[1];
                                inter_atom_list->vinter[inter_atom_list->nlocalinter][2] = atom_.v[2];
                            }
                            inter_atom_list->nlocalinter++;
                            inter_atom_list->finter.resize(inter_atom_list->nlocalinter, std::vector<double>(3));
                            inter_atom_list->rhointer.resize(inter_atom_list->nlocalinter);
                            inter_atom_list->dfinter.resize(inter_atom_list->nlocalinter);
                        } else {
                            inter_atom_list->idinter.push_back(atom_.id);
                            inter_atom_list->typeinter.push_back(atom_.type);
                            inter_atom_list->xinter.resize(inter_atom_list->nlocalinter + 1, std::vector<double>(3));
                            inter_atom_list->xinter[inter_atom_list->nlocalinter][0] = atom_.x[0];
                            inter_atom_list->xinter[inter_atom_list->nlocalinter][1] = atom_.x[1];
                            inter_atom_list->xinter[inter_atom_list->nlocalinter][2] = atom_.x[2];
                            inter_atom_list->vinter.resize(inter_atom_list->nlocalinter + 1, std::vector<double>(3));
                            inter_atom_list->vinter[inter_atom_list->nlocalinter][0] = atom_.v[0];
                            inter_atom_list->vinter[inter_atom_list->nlocalinter][1] = atom_.v[1];
                            inter_atom_list->vinter[inter_atom_list->nlocalinter][2] = atom_.v[2];
                            inter_atom_list->nlocalinter++;
                            inter_atom_list->finter.resize(inter_atom_list->nlocalinter, std::vector<double>(3));
                            inter_atom_list->rhointer.resize(inter_atom_list->nlocalinter);
                            inter_atom_list->dfinter.resize(inter_atom_list->nlocalinter);
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
    for (int i = 0; i < inter_atom_list->nlocalinter; i++) {
        if (inter_atom_list->xinter[i][0] < p_domain->getMeasuredGlobalBoxCoordLower(0)) {
            inter_atom_list->xinter[i][0] += p_domain->getMeasuredGlobalLength(0);
        } else if (inter_atom_list->xinter[i][0] >= p_domain->getMeasuredGlobalBoxCoordUpper(0)) {
            inter_atom_list->xinter[i][0] -= p_domain->getMeasuredGlobalLength(0);
        }
        if (inter_atom_list->xinter[i][1] < p_domain->getMeasuredGlobalBoxCoordLower(1)) {
            inter_atom_list->xinter[i][1] += p_domain->getMeasuredGlobalLength(1);
        } else if (inter_atom_list->xinter[i][1] >= p_domain->getMeasuredGlobalBoxCoordUpper(1)) {
            inter_atom_list->xinter[i][1] -= p_domain->getMeasuredGlobalLength(1);
        }
        if (inter_atom_list->xinter[i][2] < p_domain->getMeasuredGlobalBoxCoordLower(1)) {
            inter_atom_list->xinter[i][2] += p_domain->getMeasuredGlobalLength(2);
        } else if (inter_atom_list->xinter[i][2] >= p_domain->getMeasuredGlobalBoxCoordUpper(1)) {
            inter_atom_list->xinter[i][2] -= p_domain->getMeasuredGlobalLength(2);
        }
    }

    //判断，如果跑出晶格点的?佑峙芑鼐Ц竦悖蚍呕鼐Ц竦闶榇娲⑵湫畔?
    for (int i = 0; i < inter_atom_list->nlocalinter; i++) {
        int j, k, l;
        xtemp = inter_atom_list->xinter[i][0];
        ytemp = inter_atom_list->xinter[i][1];
        ztemp = inter_atom_list->xinter[i][2];
        j = xtemp * 2 / _latticeconst + 0.5;
        k = ytemp * 2 / _latticeconst + 0.5;
        l = ztemp * 2 / _latticeconst + 0.5;
        k = k / 2;
        l = l / 2;
        j -= p_domain->getGlobalGhostLatticeCoordLower(0);
        k -= p_domain->getGlobalGhostLatticeCoordLower(1);
        l -= p_domain->getGlobalGhostLatticeCoordLower(2);
        //判断是否在所表示晶格范围内
        if (j <= (p_domain->getSubBoxLatticeSize(0) + 2 * (ceil(_cutoffRadius / _latticeconst) + 1))
            && k <= (p_domain->getSubBoxLatticeSize(1) + (ceil(_cutoffRadius / _latticeconst) + 1))
            && l <= (p_domain->getSubBoxLatticeSize(2) + (ceil(_cutoffRadius / _latticeconst) + 1))) {
            j = atom_list->IndexOf3DIndex(j, k, l);
            AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);
            if (atom_.isInterElement()) {
                atom_.id = inter_atom_list->idinter[i];
                atom_.type = inter_atom_list->typeinter[i];
                atom_.x[0] = inter_atom_list->xinter[i][0];
                atom_.x[1] = inter_atom_list->xinter[i][1];
                atom_.x[2] = inter_atom_list->xinter[i][2];
                atom_.v[0] = inter_atom_list->vinter[i][0];
                atom_.v[1] = inter_atom_list->vinter[i][1];
                atom_.v[2] = inter_atom_list->vinter[i][2];

                inter_atom_list->idinter[i] = inter_atom_list->idinter[inter_atom_list->nlocalinter - 1];
                inter_atom_list->typeinter[i] = inter_atom_list->typeinter[inter_atom_list->nlocalinter - 1];
                inter_atom_list->xinter[i][0] = inter_atom_list->xinter[inter_atom_list->nlocalinter - 1][0];
                inter_atom_list->xinter[i][1] = inter_atom_list->xinter[inter_atom_list->nlocalinter - 1][1];
                inter_atom_list->xinter[i][2] = inter_atom_list->xinter[inter_atom_list->nlocalinter - 1][2];
                inter_atom_list->vinter[i][0] = inter_atom_list->vinter[inter_atom_list->nlocalinter - 1][0];
                inter_atom_list->vinter[i][1] = inter_atom_list->vinter[inter_atom_list->nlocalinter - 1][1];
                inter_atom_list->vinter[i][2] = inter_atom_list->vinter[inter_atom_list->nlocalinter - 1][2];

                i--;
                inter_atom_list->nlocalinter--;
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
    for (int i = 0; i < inter_atom_list->finter.size(); i++) {
        inter_atom_list->finter[i][0] = 0;
        inter_atom_list->finter[i][1] = 0;
        inter_atom_list->finter[i][2] = 0;
        inter_atom_list->rhointer[i] = 0;
    }
}

void atom::computeEam(eam *pot, Domain *domain, double &comm) {
    double starttime, stoptime;
    double xtemp, ytemp, ztemp;
    double delx, dely, delz;
    std::vector<long int>::iterator neighbourOffsetsIter;
    _type_atom_index n;
    double dist2;
    double rhoTmp, dfEmbed;
    double fpair;
    _type_atom_index kk;
    int xstart = p_domain->getGhostLatticeSize(0);
    int ystart = p_domain->getGhostLatticeSize(1);
    int zstart = p_domain->getGhostLatticeSize(2);

    // 本地晶格点上的原子计算电子云密度
    if (isAccelerateSupport()) {
//     fixme  accelerateEamRhoCalc(&(rho_spline->n), atom_list, &_cutoffRadius,
//                             &(rho_spline->invDx), rho_spline->values); // fixme
    } else { // calculate electron density use cpu only.
        for (int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
            for (int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
                for (int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                    kk = atom_list->IndexOf3DIndex(i, j, k);
                    AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                    xtemp = atom_.x[0];
                    ytemp = atom_.x[1];
                    ztemp = atom_.x[2];
                    if (!atom_.isInterElement()) {
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
                                rhoTmp = pot->rhoContribution(atom_n._tp, dist2);
                                atom_.rho += rhoTmp;
                                atom_n.rho += rhoTmp; // fixme what if the type of two atoms is different?
                            }
                        }
                    }
                }
            }
        }
    }

    //间隙原子电子云密度
    int j, k, l;
    for (int i = 0; i < inter_atom_list->nlocalinter; i++) {
        xtemp = inter_atom_list->xinter[i][0];
        ytemp = inter_atom_list->xinter[i][1];
        ztemp = inter_atom_list->xinter[i][2];
        j = xtemp * 2 / _latticeconst + 0.5;
        k = ytemp * 2 / _latticeconst + 0.5;
        l = ztemp * 2 / _latticeconst + 0.5;
        k = k / 2;
        l = l / 2;
        j -= p_domain->getGlobalGhostLatticeCoordLower(0);
        k -= p_domain->getGlobalGhostLatticeCoordLower(1);
        l -= p_domain->getGlobalGhostLatticeCoordLower(2);
        j = atom_list->IndexOf3DIndex(j, k, l);
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j);

        delx = xtemp - atom_.x[0];
        dely = ytemp - atom_.x[1];
        delz = ztemp - atom_.x[2];
        dist2 = delx * delx + dely * dely + delz * delz;
        if (dist2 < (_cutoffRadius * _cutoffRadius)) {
            rhoTmp = pot->rhoContribution(atom_._tp, dist2);
            inter_atom_list->rhointer[i] += rhoTmp;
            atom_.rho += rhoTmp; // fixme??
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
                rhoTmp = pot->rhoContribution(atom_n._tp, dist2);
                inter_atom_list->rhointer[i] += rhoTmp;
                atom_n.rho += rhoTmp; // fixme different type.
            }
            n = (j - *neighbourOffsetsIter);
            AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(n);
            delx = xtemp - atom_n.x[0];
            dely = ytemp - atom_n.x[1];
            delz = ztemp - atom_n.x[2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                rhoTmp = pot->rhoContribution(atom_._tp, dist2);
                inter_atom_list->rhointer[i] += rhoTmp;
                atom_n.rho += rhoTmp; // fixme different type.
            }
        }
        //对间隙原子遍历
        for (int k = i + 1; k < (inter_atom_list->nghostinter + inter_atom_list->nlocalinter); k++) {
            delx = xtemp - inter_atom_list->xinter[k][0];
            dely = ytemp - inter_atom_list->xinter[k][1];
            delz = ztemp - inter_atom_list->xinter[k][2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                rhoTmp = pot->rhoContribution(atom_type::Fe, dist2); // fixme, type not right
                inter_atom_list->rhointer[i] += rhoTmp;
                inter_atom_list->rhointer[k] += rhoTmp;
            }
        }
        //计算间隙原子嵌入能导数
        dfEmbed = pot->embedEnergyContribution(atom_type::Fe, inter_atom_list->rhointer[i]); // fixme atom type.
        inter_atom_list->dfinter[i] = dfEmbed;
    }

//    ofstream outfile;
    /* char tmp[20];
    sprintf(tmp, "electron_density.atom");
    outfile.open(tmp);
    for(int k =0; k < p_domain->getSubBoxLatticeSize(2) ; k++){
            for(int j = 0; j < p_domain->getSubBoxLatticeSize(1); j++){
                    for(int i =0; i < p_domain->getSubBoxLatticeSize(0) ; i++){
                             AtomElement &atom_ = atom_list->getAtomEleBySubBoxIndex(i,j,k);
                            if(!atom_.isInterElement())
                                    outfile << atom_.electron_density << std::endl;
                    }
            }
    }
for(int i = 0; i < rho_spline->n; i++){ // 1.todo remove start.
    outfile << rho_spline->spline[i][6] << std::endl;
} // 1. todo remove end.
    outfile.close();*/

    // 发送电子云密度
    // 将ghost区域的粒子的电子云密度发送给其所在的进程，得到完整的电子云密度
    starttime = MPI_Wtime();
    domain->sendrho(this);
    stoptime = MPI_Wtime();
    comm = stoptime - starttime;

    /*sprintf(tmp, "rho2.atom");
    outfile;
    outfile.open(tmp);

    for(int k =0; k < p_domain->getSubBoxLatticeSize(2) ; k++){
            for(int j = 0; j < p_domain->getSubBoxLatticeSize(1); j++){
                    for(int i =0; i < p_domain->getSubBoxLatticeSize(0) ; i++){
                             AtomElement &atom_ = atom_list->getAtomEleBySubBoxIndex(i,j,k);
                            if(!atom_.isInterElement())
                                    outfile << atom_.electron_density << std::endl;
                    }
            }
    }
    outfile.close();*/

    //本地晶格点计算嵌入能导数
    if (isAccelerateSupport()) {
//       fixme accelerateEamDfCalc(&(f_spline->n), atom_list, &_cutoffRadius,
//                            &(f_spline->invDx), f_spline->values);
    } else {
        for (int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
            for (int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
                for (int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                    kk = atom_list->IndexOf3DIndex(i, j, k);
                    AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                    dfEmbed = pot->embedEnergyContribution(atom_type::Fe, atom_.rho); // fixme atom type.
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

    // 发送嵌入能导数
    // 将本地box属于邻居进程ghost区域的粒子的嵌入能导数发送给邻居进程
    starttime = MPI_Wtime();
    domain->sendDfEmbed(this);
    stoptime = MPI_Wtime();
    comm += stoptime - starttime;

    if (isAccelerateSupport()) {
//    fixme    accelerateEamForceCalc(nullptr, atom_list, &_cutoffRadius,
//                               nullptr, nullptr, rho_spline->values);
    } else {
        /*sprintf(tmp, "f.atom");
        outfile.open(tmp);

    for(int k =0; k < p_domain->getSubBoxLatticeSize(2) ; k++){
            for(int j = 0; j < p_domain->getSubBoxLatticeSize(1); j++){
                    for(int i =0; i < p_domain->getSubBoxLatticeSize(0) ; i++){
                             AtomElement &atom_ = atom_list->getAtomEleBySubBoxIndex(i,j,k);
                                if(!atom_.isInterElement())
                                        outfile << f[kk*3] << std::endl;
                        }
                }
        }
        outfile.close();*/

        for (int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
            for (int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
                for (int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                    kk = atom_list->IndexOf3DIndex(i, j, k);
                    AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                    xtemp = atom_.x[0];
                    ytemp = atom_.x[1];
                    ztemp = atom_.x[2];
                    if (!atom_.isInterElement()) {
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
                                // fixme
                                fpair = pot->toForce(atom_._tp, atom_n._tp, dist2, atom_.df + atom_n.df);

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
    for (int i = 0; i < inter_atom_list->nlocalinter; i++) {
        int j, k, l;
        xtemp = inter_atom_list->xinter[i][0];
        ytemp = inter_atom_list->xinter[i][1];
        ztemp = inter_atom_list->xinter[i][2];
        j = xtemp * 2 / _latticeconst + 0.5;
        k = ytemp * 2 / _latticeconst + 0.5;
        l = ztemp * 2 / _latticeconst + 0.5;
        k = k / 2;
        l = l / 2;
        j -= p_domain->getGlobalGhostLatticeCoordLower(0);
        k -= p_domain->getGlobalGhostLatticeCoordLower(1);
        l -= p_domain->getGlobalGhostLatticeCoordLower(2);
        j = atom_list->IndexOf3DIndex(j, k, l);
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(j); // cgs: 间隙原子所在晶格处的原子

        delx = xtemp - atom_.x[0];
        dely = ytemp - atom_.x[1];
        delz = ztemp - atom_.x[2];
        dist2 = delx * delx + dely * dely + delz * delz;
        if (dist2 < (_cutoffRadius * _cutoffRadius)) {
            // fixme inter atom type
            fpair = pot->toForce(atom_type::Fe, atom_._tp, dist2, inter_atom_list->dfinter[i] + atom_.df);

            inter_atom_list->finter[i][0] += delx * fpair;
            inter_atom_list->finter[i][1] += dely * fpair;
            inter_atom_list->finter[i][2] += delz * fpair;

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
                // fixme
                fpair = pot->toForce(atom_type::Fe, atom_n._tp, dist2, inter_atom_list->dfinter[i] + atom_n.df);

                inter_atom_list->finter[i][0] += delx * fpair;
                inter_atom_list->finter[i][1] += dely * fpair;
                inter_atom_list->finter[i][2] += delz * fpair;

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
                // fixme
                fpair = pot->toForce(atom_type::Fe, atom_n._tp, dist2, inter_atom_list->dfinter[i] + atom_.df);

                inter_atom_list->finter[i][0] += delx * fpair;
                inter_atom_list->finter[i][1] += dely * fpair;
                inter_atom_list->finter[i][2] += delz * fpair;

                atom_.f[0] -= delx * fpair;
                atom_.f[1] -= dely * fpair;
                atom_.f[2] -= delz * fpair;
            }
        }
        //对间隙原子遍历
        for (int k = i + 1; k < (inter_atom_list->nghostinter + inter_atom_list->nlocalinter); k++) {
            delx = xtemp - inter_atom_list->xinter[k][0];
            dely = ytemp - inter_atom_list->xinter[k][1];
            delz = ztemp - inter_atom_list->xinter[k][2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                // fixme
                fpair = pot->toForce(atom_type::Fe, atom_type::Fe, dist2,
                                     inter_atom_list->dfinter[i] + inter_atom_list->dfinter[k]);

                inter_atom_list->finter[i][0] += delx * fpair;
                inter_atom_list->finter[i][1] += dely * fpair;
                inter_atom_list->finter[i][2] += delz * fpair;

                inter_atom_list->finter[k][0] -= delx * fpair;
                inter_atom_list->finter[k][1] -= dely * fpair;
                inter_atom_list->finter[k][2] -= delz * fpair;
            }
        }
    }
}

unsigned long atom::getinteridsendsize() {
    return interbuf.size();
}

void atom::getatomx(int direction, std::vector<std::vector<_type_atom_id> > &sendlist) {
    _type_atom_id i;
    if (direction == 0) {
        //找到要发送到邻居进程的区域
        int xstart = p_domain->getGhostLatticeSize(0);
        int ystart = p_domain->getGhostLatticeSize(1);
        int zstart = p_domain->getGhostLatticeSize(2);
        int xstop = xstart + p_domain->getGhostLatticeSize(0); // note: this is ghost lattice size.
        int ystop = ystart + p_domain->getSubBoxLatticeSize(1);
        int zstop = zstart + p_domain->getSubBoxLatticeSize(2);

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = atom_list->IndexOf3DIndex(ix, iy, iz);
                    sendlist[0].push_back(i);
                }
            }
        }
    } else {
        //找到要发送到邻居进程的区域
        int xstart = p_domain->getGhostLatticeSize(0) + p_domain->getSubBoxLatticeSize(0) - ((_cutlattice) * 2);
        int ystart = p_domain->getGhostLatticeSize(1);
        int zstart = p_domain->getGhostLatticeSize(2);
        int xstop = p_domain->getGhostLatticeSize(0) + p_domain->getSubBoxLatticeSize(0);
        int ystop = ystart + p_domain->getSubBoxLatticeSize(1);
        int zstop = zstart + p_domain->getSubBoxLatticeSize(2);

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = atom_list->IndexOf3DIndex(ix, iy, iz);
                    sendlist[1].push_back(i);
                }
            }
        }
    }
}

void atom::getatomy(int direction, std::vector<std::vector<_type_atom_id> > &sendlist) {
    int i;
    if (direction == 0) {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = p_domain->getGhostLatticeSize(1);
        int zstart = p_domain->getGhostLatticeSize(2);
        int xstop = p_domain->getGhostExtLatticeSize(0);
        int ystop = ystart + _cutlattice;
        int zstop = zstart + p_domain->getSubBoxLatticeSize(2);

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = atom_list->IndexOf3DIndex(ix, iy, iz);
                    sendlist[2].push_back(i);
                }
            }
        }
    } else {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = p_domain->getGhostLatticeSize(1) + p_domain->getSubBoxLatticeSize(1) - (_cutlattice);
        int zstart = p_domain->getGhostLatticeSize(2);
        int xstop = p_domain->getGhostExtLatticeSize(0);
        int ystop = p_domain->getGhostLatticeSize(1) + p_domain->getSubBoxLatticeSize(1);
        int zstop = zstart + p_domain->getSubBoxLatticeSize(2);

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = atom_list->IndexOf3DIndex(ix, iy, iz);
                    sendlist[3].push_back(i);
                }
            }
        }
    }
}

void atom::getatomz(int direction, std::vector<std::vector<_type_atom_id> > &sendlist) {
    int i;
    if (direction == 0) {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = 0;
        int zstart = p_domain->getGhostLatticeSize(2);
        int xstop = p_domain->getGhostExtLatticeSize(0);
        int ystop = p_domain->getGhostExtLatticeSize(1);
        int zstop = zstart + _cutlattice;

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = atom_list->IndexOf3DIndex(ix, iy, iz);
                    sendlist[4].push_back(i);
                }
            }
        }
    } else {
        //找到要发送到邻居进程的区域
        int xstart = 0;
        int ystart = 0;
        int zstart = p_domain->getGhostLatticeSize(2) + p_domain->getSubBoxLatticeSize(2) - (_cutlattice);
        int xstop = p_domain->getGhostExtLatticeSize(0);
        int ystop = p_domain->getGhostExtLatticeSize(1);
        int zstop = p_domain->getGhostLatticeSize(2) + p_domain->getSubBoxLatticeSize(2);

        //要发送要邻居进程区域内的分子指针
        for (int iz = zstart; iz < zstop; iz++) {
            for (int iy = ystart; iy < ystop; iy++) {
                for (int ix = xstart; ix < xstop; ix++) {
                    i = atom_list->IndexOf3DIndex(ix, iy, iz);
                    sendlist[5].push_back(i);
                }
            }
        }
    }
}

void atom::getIntertosend(int d, int direction, double ghostlengh, std::vector<int> &sendlist) {
    double low, high;
    if (d == 0) {
        if (direction == 0) {
            low = p_domain->getMeasuredSubBoxLowerBounding(0);
            high = p_domain->getMeasuredSubBoxLowerBounding(0) + ghostlengh;
            for (int i = 0; i < inter_atom_list->nlocalinter; i++) {
                if (inter_atom_list->xinter[i][0] < high && inter_atom_list->xinter[i][0] >= low) {
                    sendlist.push_back(i);
                }
            }
        } else {
            low = p_domain->getMeasuredSubBoxUpperBounding(0) - ghostlengh;
            high = p_domain->getMeasuredSubBoxUpperBounding(0);
            for (int i = 0; i < inter_atom_list->nlocalinter; i++) {
                if (inter_atom_list->xinter[i][0] <= high && inter_atom_list->xinter[i][0] > low) {
                    sendlist.push_back(i);
                }
            }
        }
    } else if (d == 1) {
        if (direction == 0) {
            low = p_domain->getMeasuredSubBoxLowerBounding(1);
            high = p_domain->getMeasuredSubBoxLowerBounding(1) + ghostlengh;
            for (int i = 0; i < inter_atom_list->nlocalinter + inter_atom_list->nghostinter; i++) {
                if (inter_atom_list->xinter[i][1] < high && inter_atom_list->xinter[i][1] >= low) {
                    sendlist.push_back(i);
                }
            }
        } else {
            low = p_domain->getMeasuredSubBoxUpperBounding(1) - ghostlengh;
            high = p_domain->getMeasuredSubBoxUpperBounding(1);
            for (int i = 0; i < inter_atom_list->nlocalinter + inter_atom_list->nghostinter; i++) {
                if (inter_atom_list->xinter[i][1] <= high && inter_atom_list->xinter[i][1] > low) {
                    sendlist.push_back(i);
                }
            }
        }
    } else {
        if (direction == 0) {
            low = p_domain->getMeasuredSubBoxLowerBounding(2);
            high = p_domain->getMeasuredSubBoxLowerBounding(2) + ghostlengh;
            for (int i = 0; i < inter_atom_list->nlocalinter + inter_atom_list->nghostinter; i++) {
                if (inter_atom_list->xinter[i][2] < high && inter_atom_list->xinter[i][2] >= low) {
                    sendlist.push_back(i);
                }
            }
        } else {
            low = p_domain->getMeasuredSubBoxUpperBounding(2) - ghostlengh;
            high = p_domain->getMeasuredSubBoxUpperBounding(2);
            for (int i = 0; i < inter_atom_list->nlocalinter + inter_atom_list->nghostinter; i++) {
                if (inter_atom_list->xinter[i][2] <= high && inter_atom_list->xinter[i][2] > low) {
                    sendlist.push_back(i);
                }
            }
        }
    }
}

int atom::getintersendnum(int dimension, int direction) {
    interbuf.clear();
    for (int i = 0; i < inter_atom_list->nlocalinter; i++) {
        if (direction == 0) {
            // we assume that, a atom cannot cross 2 or more than 2 sub-boxes
            if (inter_atom_list->xinter[i][dimension] < p_domain->getMeasuredSubBoxLowerBounding(dimension)) {
                interbuf.push_back(i);
            }
        } else {
            if (inter_atom_list->xinter[i][dimension] >= p_domain->getMeasuredSubBoxUpperBounding(dimension)) {
                interbuf.push_back(i);
            }
        }
    }
    return interbuf.size();
}

void atom::computefirst(double dtInv2m, double dt) {
    //本地晶格点上的原子求解运动方程第一步
    atom_list->foreachSubBoxAtom(
            [dtInv2m, dt](AtomElement &_atom_ref) {
                if (!_atom_ref.isInterElement()) {
                    _atom_ref.v[0] = _atom_ref.v[0] + dtInv2m * _atom_ref.f[0];
                    _atom_ref.v[1] = _atom_ref.v[1] + dtInv2m * _atom_ref.f[1];
                    _atom_ref.v[2] = _atom_ref.v[2] + dtInv2m * _atom_ref.f[2];
                    _atom_ref.x[0] += dt * _atom_ref.v[0];
                    _atom_ref.x[1] += dt * _atom_ref.v[1];
                    _atom_ref.x[2] += dt * _atom_ref.v[2];
                }
            }
    );

    //本地间隙原子求解运动方程第一步
    for (int i = 0; i < inter_atom_list->nlocalinter; i++) {
        for (unsigned short d = 0; d < 3; ++d) {
            inter_atom_list->vinter[i][d] = inter_atom_list->vinter[i][d] + dtInv2m * inter_atom_list->finter[i][d];
            inter_atom_list->xinter[i][d] += dt * inter_atom_list->vinter[i][d];
        }
    }
}

void atom::computesecond(double dtInv2m) {
    //本地晶格点上的原子求解运动方程第二步
    atom_list->foreachSubBoxAtom(
            [dtInv2m](AtomElement &_atom_ref) {
                // fixme, excluding InterElement. (add if: isInterElement)
                for (unsigned short d = 0; d < DIMENSION; ++d) {
                    _atom_ref.v[d] += dtInv2m * _atom_ref.f[d];
                }
            }
    );
    //本地间隙原子求解运动方程第二步
    for (int i = 0; i < inter_atom_list->nlocalinter; i++) {
        for (unsigned short d = 0; d < 3; ++d) {
            inter_atom_list->vinter[i][d] = inter_atom_list->vinter[i][d] + dtInv2m * inter_atom_list->finter[i][d];
        }
    }
}

void atom::print_force() {
    char tmp[20];
    sprintf(tmp, "force.txt");
    std::ofstream outfile;
    outfile.open(tmp);

    atom_list->foreachSubBoxAtom(
            [&outfile](AtomElement &_atom_ref) {
                outfile << _atom_ref.f[0] << " " << _atom_ref.f[1] << " " << _atom_ref.f[2] << std::endl;
            }
    );

    long kk;
    int xstart = p_domain->getGhostLatticeSize(0);
    int ystart = p_domain->getGhostLatticeSize(1);
    int zstart = p_domain->getGhostLatticeSize(2);
    std::cout << "print_force" << std::endl;
    for (int k = zstart; k < p_domain->getSubBoxLatticeSize(2) + zstart; k++) {
        for (int j = ystart; j < p_domain->getSubBoxLatticeSize(1) + ystart; j++) {
            for (int i = xstart; i < p_domain->getSubBoxLatticeSize(0) + xstart; i++) {
                kk = atom_list->IndexOf3DIndex(i, j, k);
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                outfile << atom_.f[0] << " " << atom_.f[1] << " " << atom_.f[2] << std::endl;
            }
        }
    }
    outfile.close();
}

void atom::setv(int lat[4], double collision_v[3]) {
    long kk;
    if ((lat[0] * 2) >= p_domain->getGlobalSubBoxLatticeCoordLower(0) &&
        (lat[0] * 2) < (p_domain->getGlobalSubBoxLatticeCoordLower(0) + p_domain->getSubBoxLatticeSize(0))
        && lat[1] >= p_domain->getGlobalSubBoxLatticeCoordLower(1) &&
        lat[1] < (p_domain->getGlobalSubBoxLatticeCoordLower(1) + p_domain->getSubBoxLatticeSize(1))
        && lat[2] >= p_domain->getGlobalSubBoxLatticeCoordLower(2) &&
        lat[2] < (p_domain->getGlobalSubBoxLatticeCoordLower(2) + p_domain->getSubBoxLatticeSize(2))) {
        kk = (atom_list->IndexOf3DIndex(lat[0] * 2 - p_domain->getGlobalGhostLatticeCoordLower(0),
                                        lat[1] - p_domain->getGlobalGhostLatticeCoordLower(1),
                                        lat[2] - p_domain->getGlobalGhostLatticeCoordLower(2)) + lat[4]);
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
