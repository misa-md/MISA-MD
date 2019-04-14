#include <iostream>
#include <iterator>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <utils/mpi_domain.h>
#include <logs/logs.h>
#include <eam.h>
#include <comm.hpp>

#include "atom.h"
#include "atom/ws_utils.h"
#include "pack/rho_packer.h"
#include "pack/force_packer.h"
#include "pack/df_embed_packer.h"
#include "hardware_accelerate.hpp" // use hardware(eg.GPU, MIC,Sunway slave cores.) to achieve calculate accelerating.

atom::atom(comm::Domain *domain)
        : AtomSet(domain->lattice_const * domain->cutoff_radius_factor,
                  domain->lattice_size_ghost_extended,
                  domain->lattice_size_sub_box,
                  domain->lattice_size_ghost),
          p_domain(domain) {
    if (isAccelerateSupport()) {
        accelerateInit(p_domain->dbx_lattice_coord_sub_box_region.x_low,
                       p_domain->dbx_lattice_coord_sub_box_region.y_low,
                       p_domain->dbx_lattice_coord_sub_box_region.z_low,
                       p_domain->dbx_lattice_size_sub_box[0],
                       p_domain->dbx_lattice_size_sub_box[1],
                       p_domain->dbx_lattice_size_sub_box[2],
                       p_domain->dbx_lattice_coord_ghost_region.x_low,
                       p_domain->dbx_lattice_coord_ghost_region.y_low,
                       p_domain->dbx_lattice_coord_ghost_region.z_low,
                       p_domain->dbx_lattice_size_ghost_extended[0],
                       p_domain->dbx_lattice_size_ghost_extended[1],
                       p_domain->dbx_lattice_size_ghost_extended[2]);
    }
}

int atom::decide() {
    inter_atom_list->clearGhost();
    int nflag = 0;
    long kk = 0;
    double dist;
    double xtemp, ytemp, ztemp;

    //对本地晶格点原子进行判断，看是否运动为间隙原子
    for (_type_atom_index k = 0; k < p_domain->dbx_lattice_size_sub_box[2]; k++) {
        for (_type_atom_index j = 0; j < p_domain->dbx_lattice_size_sub_box[1]; j++) {
            for (_type_atom_index i = 0; i < p_domain->dbx_lattice_size_sub_box[0]; i++) {
//                kk = atom_list->IndexOf3DIndex(i, j, k);
                AtomElement &atom_ = atom_list->getAtomEleBySubBoxIndex(i, j, k); // todo long type
                if (!atom_.isInterElement()) {
                    xtemp = (i + p_domain->dbx_lattice_coord_sub_box_region.x_low) * 0.5 * p_domain->lattice_const;
                    ytemp = (j + p_domain->dbx_lattice_coord_sub_box_region.y_low + (i % 2) * 0.5) *
                            p_domain->lattice_const;
                    ztemp = (k + p_domain->dbx_lattice_coord_sub_box_region.z_low + (i % 2) * 0.5) *
                            p_domain->lattice_const;
                    dist = (atom_.x[0] - xtemp) * (atom_.x[0] - xtemp);
                    dist += (atom_.x[1] - ytemp) * (atom_.x[1] - ytemp);
                    dist += (atom_.x[2] - ztemp) * (atom_.x[2] - ztemp);
                    if (dist > (pow(0.2 * p_domain->lattice_const, 2.0))) { /**超过距离则判断为间隙原子*/
                        inter_atom_list->addInterAtom(atom_);
                        atom_.type = atom_type::INVALID;
                        atom_.v[0] = 0;
                        atom_.v[1] = 0;
                        atom_.v[2] = 0;
                        nflag = 1;
                    }
                }
            }
        }
    }

    // 如果间隙原子跑入晶格点,且晶格点为空位, 则空位-间隙发生复合.
    // note: the out-of-lattice atoms can be out of sub-box.
    _type_inter_list::iterator inter_it;
    for (inter_it = inter_atom_list->inter_list.begin(); inter_it != inter_atom_list->inter_list.end();) {
        AtomElement &inter_ref = *inter_it;
        // get the near atom of the inter atom
        // If we use func findNearLatAtom to find a near atom of an inter atom in atoms list
        // the near atom can be an ghost atom (but the position of that ghost atom may still be in sub-box).
        // we should find near atom only in lattice atoms(exclude ghost atoms), so we use func finNearLatAtomInSubBox.
        AtomElement *near_atom = ws::findNearLatAtomInSubBox(atom_list, inter_ref, p_domain);
        // the near atom must be in sub-box, and it is in the lattice atom lists.
        if (near_atom != nullptr && near_atom->isInterElement() &&
            ws::isOutBox(*near_atom, p_domain) == box::IN_BOX) {
            near_atom->id = inter_ref.id;
            near_atom->type = inter_ref.type; // set type to valid.
            near_atom->x[0] = inter_ref.x[0];
            near_atom->x[1] = inter_ref.x[1];
            near_atom->x[2] = inter_ref.x[2];
            near_atom->v[0] = inter_ref.v[0];
            near_atom->v[1] = inter_ref.v[1];
            near_atom->v[2] = inter_ref.v[2];

            // remove this atom from inter list.
            inter_it = inter_atom_list->removeInter(inter_it);
        } else {
            inter_it++;
        }
    }
    return nflag;
}

void atom::clearForce() {
    for (_type_atom_index i = 0; i < atom_list->size(); i++) {
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(i);
        atom_.f[0] = 0;
        atom_.f[1] = 0;
        atom_.f[2] = 0;
        atom_.rho = 0;
    }
    for (AtomElement &inter_ref :inter_atom_list->inter_list) {
        inter_ref.f[0] = 0;
        inter_ref.f[1] = 0;
        inter_ref.f[2] = 0;
        inter_ref.rho = 0;
    }
}

void atom::computeEam(eam *pot, comm::Domain *domain, double &comm) {
    double starttime, stoptime;
    inter_atom_list->makeIndex(atom_list, p_domain); // create index for inter atom and inter ghost atoms.

    latRho(pot, domain, comm);
    interRho(pot, domain, comm);
//    ofstream outfile;
    /* char tmp[20];
    sprintf(tmp, "electron_density.atom");
    outfile.open(tmp);
    int j, k, l;
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

    {
        // 发送电子云密度
        // 将ghost区域的粒子的电子云密度发送给其所在的进程，得到完整的电子云密度
        starttime = MPI_Wtime();
        RhoPacker rho_packer(getAtomListRef(), atom_list->sendlist, atom_list->recvlist);
        comm::neiSendReceive<double, MPI_DOUBLE>(&rho_packer, MPIDomain::toCommProcess(), p_domain->rank_id_neighbours);
        stoptime = MPI_Wtime();
        comm = stoptime - starttime;
    }

    /*sprintf(tmp, "rho2.atom");
    outfile;
    outfile.open(tmp);
    int j, k, l;
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
    latDf(pot, domain, comm);

    /*sprintf(tmp, "df.atom");
    outfile.open(tmp);
    for(int i = 0; i < f_spline->n; i++){
        outfile << i << " " << f_spline->spline[i][6] << std::endl;
    }
    outfile.close();*/

    {
        // 发送嵌入能导数
        // 将本地box属于邻居进程ghost区域的粒子的嵌入能导数发送给邻居进程
        starttime = MPI_Wtime();
        DfEmbedPacker packer(getAtomListRef(), *inter_atom_list,
                             atom_list->sendlist, atom_list->recvlist,
                             inter_atom_list->intersendlist,
                             inter_atom_list->interrecvlist);
        comm::neiSendReceive<double, MPI_DOUBLE>(&packer, MPIDomain::toCommProcess(), p_domain->rank_id_neighbours);
        stoptime = MPI_Wtime();
        comm += stoptime - starttime;
    }

    // force for local lattice.
    latForce(pot, domain, comm);

    /*sprintf(tmp, "f.atom");  // 2.todo remove start.
      outfile.open(tmp);
      for(int i = 0; i < phi_spline->n; i++){
         outfile << i << " " << phi_spline->spline[i][6] << std::endl;
      }
      outfile.close();*/ // 2.todo remove end.
    //间隙原子计算嵌入能和对势带来的力
    interForce(pot, domain, comm);
}

void atom::latRho(eam *pot, comm::Domain *domain, double &comm) {
    double xtemp, ytemp, ztemp;
    double delx, dely, delz;
    _type_atom_index kk;
    _type_atom_index n;
    double dist2;
    int xstart = p_domain->dbx_lattice_size_ghost[0];
    int ystart = p_domain->dbx_lattice_size_ghost[1];
    int zstart = p_domain->dbx_lattice_size_ghost[2];

    // 本地晶格点上的原子计算电子云密度
    if (isAccelerateSupport()) {
//     fixme  accelerateEamRhoCalc(&(rho_spline->n), atom_list, &_cutoffRadius,
//                             &(rho_spline->invDx), rho_spline->values); // fixme
    } else { // calculate electron density use cpu only.
        for (int k = zstart; k < p_domain->dbx_lattice_size_sub_box[2] + zstart; k++) {
            for (int j = ystart; j < p_domain->dbx_lattice_size_sub_box[1] + ystart; j++) {
                for (int i = xstart; i < p_domain->dbx_lattice_size_sub_box[0] + xstart; i++) {
                    kk = atom_list->IndexOf3DIndex(i, j, k);
                    AtomElement &atom_central = atom_list->getAtomEleByLinearIndex(kk);
                    xtemp = atom_central.x[0];
                    ytemp = atom_central.x[1];
                    ztemp = atom_central.x[2];
                    if (!atom_central.isInterElement()) {
                        //对晶格点邻居原子遍历
                        // only consider the atoms whose id is bigger than {@var atom_central}, just single side.
                        AtomNei::iterator nei_itl_end = neighbours->end(true, i, j, k);
                        for (AtomNei::iterator nei_itl = neighbours->begin(true, i, j, k);
                             nei_itl != nei_itl_end; ++nei_itl) {
                            AtomElement &atom_neighbour = *nei_itl;
                            if (atom_neighbour.isInterElement()) {
                                continue;
                            }
                            delx = xtemp - atom_neighbour.x[0];
                            dely = ytemp - atom_neighbour.x[1];
                            delz = ztemp - atom_neighbour.x[2];
                            dist2 = delx * delx + dely * dely + delz * delz;
                            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                                atom_central.rho += pot->rhoContribution(
                                        atom_type::getTypeIdByType(atom_neighbour.type), dist2);
                                atom_neighbour.rho += pot->rhoContribution(
                                        atom_type::getTypeIdByType(atom_central.type), dist2);
                                // fixme
                            }
                        }
                    }
                }
            }
        }
    }
}

void atom::interRho(eam *pot, comm::Domain *domain, double &comm) {
    double delx, dely, delz;
    _type_atom_index n;
    double dist2;
    double dfEmbed;

    std::vector<_type_atom_index>::iterator neighbourOffsetsIter;

    //间隙原子电子云密度
    _type_atom_index near_atom_index;
    // note that some inter atoms may be not in sub-box, we should exclude them.
    for (_type_inter_list::iterator inter_it = inter_atom_list->inter_list.begin();
         inter_it != inter_atom_list->inter_list.end(); inter_it++) {
        // get index of nearest atom of inter atoms.
        near_atom_index = ws::findNearLatIndexInSubBox(atom_list, *inter_it, p_domain);
        if (near_atom_index == box::IndexNotExists) {
            continue; // make sure the inter atoms is in sub box.
            // todo find a good way to filter out-of-box atoms while exchanging inter atoms.
        }
        AtomElement &atom_near = atom_list->getAtomEleByLinearIndex(near_atom_index);
        delx = (*inter_it).x[0] - atom_near.x[0];
        dely = (*inter_it).x[1] - atom_near.x[1];
        delz = (*inter_it).x[2] - atom_near.x[2];
        dist2 = delx * delx + dely * dely + delz * delz;
        if (!atom_near.isInterElement() && dist2 < (_cutoffRadius * _cutoffRadius)) {
            inter_it->rho += pot->rhoContribution(atom_type::getTypeIdByType(atom_near.type), dist2);
            atom_near.rho += pot->rhoContribution(atom_type::getTypeIdByType((*inter_it).type), dist2);
            // fixme
        }

        _type_atom_index x, y, z;
        atom_list->get3DIndexByLinearIndex(near_atom_index, x, y, z);
        // rho between inter atoms and lattice atoms (use full neighbour index).
        AtomNei::iterator nei_full_itl_end = neighbours->end(false, x, y, z);
        for (AtomNei::iterator nei_itl = neighbours->begin(false, x, y, z);
             nei_itl != nei_full_itl_end; ++nei_itl) {
            AtomElement &lat_nei_atom = *nei_itl; // this is a lattice atom.
            if (!lat_nei_atom.isInterElement()) {
                delx = (*inter_it).x[0] - lat_nei_atom.x[0];
                dely = (*inter_it).x[1] - lat_nei_atom.x[1];
                delz = (*inter_it).x[2] - lat_nei_atom.x[2];
                dist2 = delx * delx + dely * dely + delz * delz;
                if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                    (*inter_it).rho += pot->rhoContribution(
                            atom_type::getTypeIdByType(lat_nei_atom.type), dist2);
                    lat_nei_atom.rho += pot->rhoContribution(
                            atom_type::getTypeIdByType((*inter_it).type), dist2);
                    // fixme
                }
            }
                }

        // rho between inter atoms and inter atoms (use half neighbour index).
        AtomNei::iterator nei_half_itl_end = neighbours->end(true, x, y, z);
        for (AtomNei::iterator nei_itl = neighbours->begin(true, x, y, z);
             nei_itl != nei_half_itl_end; ++nei_itl) {
            const _type_atom_index inter_nei_id = atom_list->IndexOf3DIndex(
                    nei_itl.cur_index_x, nei_itl.cur_index_y,
                    nei_itl.cur_index_z); // get index of the neighbour lattice.
            // get intel atoms on this neighbour lattice and calculate inter-rho.
            inter_map_range inter_map_range = inter_atom_list->inter_map.equal_range(inter_nei_id);
            for (inter_map_range_itl itl = inter_map_range.first; itl != inter_map_range.second; ++itl) {
                delx = (*inter_it).x[0] - itl->second->x[0];
                dely = (*inter_it).x[1] - itl->second->x[1];
                delz = (*inter_it).x[2] - itl->second->x[2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                    (*inter_it).rho += pot->rhoContribution(atom_type::getTypeIdByType(itl->second->type), dist2);
                    itl->second->rho += pot->rhoContribution(atom_type::getTypeIdByType((*inter_it).type), dist2);
                }
            }
        }
        // todo inter ghost atoms -> cell atoms
        //计算间隙原子嵌入能导数
        // fixme
        dfEmbed = pot->embedEnergyContribution(atom_type::getTypeIdByType((*inter_it).type), (*inter_it).rho);
        (*inter_it).df = dfEmbed;
    }
}

void atom::latDf(eam *pot, comm::Domain *domain, double &comm) {
    double dfEmbed;
    int xstart = p_domain->dbx_lattice_size_ghost[0];
    int ystart = p_domain->dbx_lattice_size_ghost[1];
    int zstart = p_domain->dbx_lattice_size_ghost[2];
    _type_atom_index kk;

    //本地晶格点计算嵌入能导数
    if (isAccelerateSupport()) {
//       fixme accelerateEamDfCalc(&(f_spline->n), atom_list, &_cutoffRadius,
//                            &(f_spline->invDx), f_spline->values);
    } else {
        for (int k = zstart; k < p_domain->dbx_lattice_size_sub_box[2] + zstart; k++) {
            for (int j = ystart; j < p_domain->dbx_lattice_size_sub_box[1] + ystart; j++) {
                for (int i = xstart; i < p_domain->dbx_lattice_size_sub_box[0] + xstart; i++) {
                    kk = atom_list->IndexOf3DIndex(i, j, k);
                    AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                    if (atom_.isInterElement()) {
                        continue;
                    }
                    dfEmbed = pot->embedEnergyContribution(atom_type::getTypeIdByType(atom_.type), atom_.rho); // fixme
                    atom_.df = dfEmbed;
                }
            }
        }
    }
}

void atom::latForce(eam *pot, comm::Domain *domain, double &comm) {
    double xtemp, ytemp, ztemp;
    double delx, dely, delz;
    _type_atom_index kk;
    _type_atom_index n;
    double dist2;
    double fpair;

    int xstart = p_domain->dbx_lattice_size_ghost[0];
    int ystart = p_domain->dbx_lattice_size_ghost[1];
    int zstart = p_domain->dbx_lattice_size_ghost[2];

    std::vector<_type_atom_index>::iterator neighbourOffsetsIter;

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

        for (int k = zstart; k < p_domain->dbx_lattice_size_sub_box[2] + zstart; k++) {
            for (int j = ystart; j < p_domain->dbx_lattice_size_sub_box[1] + ystart; j++) {
                for (int i = xstart; i < p_domain->dbx_lattice_size_sub_box[0] + xstart; i++) {
                    kk = atom_list->IndexOf3DIndex(i, j, k);
                    AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                    xtemp = atom_.x[0];
                    ytemp = atom_.x[1];
                    ztemp = atom_.x[2];
                    if (atom_.isInterElement()) {
                        continue;
                    }
                    //对晶格点邻居原子遍历
                    for (neighbourOffsetsIter = NeighbourOffsets.begin();
                         neighbourOffsetsIter != NeighbourOffsets.end(); neighbourOffsetsIter++) {
                        n = (kk + *neighbourOffsetsIter); // todo what it is inter atom?
                        AtomElement &atom_n = atom_list->getAtomEleByLinearIndex(n);
                        delx = xtemp - atom_n.x[0];
                        dely = ytemp - atom_n.x[1];
                        delz = ztemp - atom_n.x[2];
                        dist2 = delx * delx + dely * dely + delz * delz;
                        if (dist2 < (_cutoffRadius * _cutoffRadius) && !atom_n.isInterElement()) {
                            // fixme
                            fpair = pot->toForce(atom_type::getTypeIdByType(atom_.type),
                                                 atom_type::getTypeIdByType(atom_n.type),
                                                 dist2, atom_.df + atom_n.df);

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
    } // end of if-isAccelerateSupport.
}

void atom::interForce(eam *pot, comm::Domain *domain, double &comm) {
    double delx, dely, delz;
    _type_atom_index n;
    double dist2;
    double fpair;

    std::vector<_type_atom_index>::iterator neighbourOffsetsIter;

    //间隙原子计算嵌入能和对势带来的力
    _type_atom_index _atom_near_index;
    for (_type_inter_list::iterator inter_it = inter_atom_list->inter_list.begin();
         inter_it != inter_atom_list->inter_list.end(); inter_it++) {
        _atom_near_index = ws::findNearLatIndexInSubBox(atom_list, *inter_it, p_domain);
        if (_atom_near_index == box::IndexNotExists) {
            continue; // make sure the inter atoms is in sub box.
            // todo find a good way to filter out-of-box atoms while exchanging inter atoms.
        }
        // 间隙原子所在晶格处的原子
        AtomElement &atom_central = atom_list->getAtomEleByLinearIndex(_atom_near_index);

        delx = (*inter_it).x[0] - atom_central.x[0];
        dely = (*inter_it).x[1] - atom_central.x[1];
        delz = (*inter_it).x[2] - atom_central.x[2];
        dist2 = delx * delx + dely * dely + delz * delz;
        if (dist2 < (_cutoffRadius * _cutoffRadius) && !atom_central.isInterElement()) {
            // fixme
            fpair = pot->toForce(
                    atom_type::getTypeIdByType((*inter_it).type),
                    atom_type::getTypeIdByType(atom_central.type),
                    dist2, (*inter_it).df + atom_central.df);

            (*inter_it).f[0] += delx * fpair;
            (*inter_it).f[1] += dely * fpair;
            (*inter_it).f[2] += delz * fpair;

            atom_central.f[0] -= delx * fpair;
            atom_central.f[1] -= dely * fpair;
            atom_central.f[2] -= delz * fpair;
        }
        for (neighbourOffsetsIter = NeighbourOffsets.begin();
             neighbourOffsetsIter != NeighbourOffsets.end(); neighbourOffsetsIter++) {
            //计算间隙原子的所有邻居
            // up lattice atoms
            n = (_atom_near_index + *neighbourOffsetsIter);
            AtomElement &atom_neighbour_up = atom_list->getAtomEleByLinearIndex(n);
            delx = (*inter_it).x[0] - atom_neighbour_up.x[0];
            dely = (*inter_it).x[1] - atom_neighbour_up.x[1];
            delz = (*inter_it).x[2] - atom_neighbour_up.x[2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius) && !atom_neighbour_up.isInterElement()) {
                // fixme
                fpair = pot->toForce(
                        atom_type::getTypeIdByType((*inter_it).type),
                        atom_type::getTypeIdByType(atom_neighbour_up.type),
                        dist2, (*inter_it).df + atom_neighbour_up.df);

                (*inter_it).f[0] += delx * fpair;
                (*inter_it).f[1] += dely * fpair;
                (*inter_it).f[2] += delz * fpair;

                atom_neighbour_up.f[0] -= delx * fpair;
                atom_neighbour_up.f[1] -= dely * fpair;
                atom_neighbour_up.f[2] -= delz * fpair;
            }
            // up inter atoms
            {
                inter_map_range inter_map_range_up = inter_atom_list->inter_map.equal_range(n);
                for (inter_map_range_itl itl_up = inter_map_range_up.first;
                     itl_up != inter_map_range_up.second; ++itl_up) {
                    delx = (*inter_it).x[0] - itl_up->second->x[0];
                    dely = (*inter_it).x[1] - itl_up->second->x[1];
                    delz = (*inter_it).x[2] - itl_up->second->x[2];
                    dist2 = delx * delx + dely * dely + delz * delz;
                    if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                        fpair = pot->toForce(
                                atom_type::getTypeIdByType((*inter_it).type),
                                atom_type::getTypeIdByType(itl_up->second->type),
                                dist2, (*inter_it).df + itl_up->second->df);

                        (*inter_it).f[0] += delx * fpair;
                        (*inter_it).f[1] += dely * fpair;
                        (*inter_it).f[2] += delz * fpair;

                        itl_up->second->f[0] -= delx * fpair;
                        itl_up->second->f[1] -= dely * fpair;
                        itl_up->second->f[2] -= delz * fpair;
                    }
                }
            }
            // down lattice atoms
            n = (_atom_near_index - *neighbourOffsetsIter);
            AtomElement &atom_neighbour_down = atom_list->getAtomEleByLinearIndex(n);

            delx = (*inter_it).x[0] - atom_neighbour_down.x[0];
            dely = (*inter_it).x[1] - atom_neighbour_down.x[1];
            delz = (*inter_it).x[2] - atom_neighbour_down.x[2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius) && !atom_neighbour_down.isInterElement()) {
                // fixme
                fpair = pot->toForce(
                        atom_type::getTypeIdByType((*inter_it).type),
                        atom_type::getTypeIdByType(atom_neighbour_down.type),
                        dist2, (*inter_it).df + atom_neighbour_down.df);

                (*inter_it).f[0] += delx * fpair;
                (*inter_it).f[1] += dely * fpair;
                (*inter_it).f[2] += delz * fpair;

                atom_neighbour_down.f[0] -= delx * fpair;
                atom_neighbour_down.f[1] -= dely * fpair;
                atom_neighbour_down.f[2] -= delz * fpair;
            }
            // down inter atoms
            {
                inter_map_range inter_map_range_up = inter_atom_list->inter_map.equal_range(n);
                for (inter_map_range_itl itl_up = inter_map_range_up.first;
                     itl_up != inter_map_range_up.second; ++itl_up) {
                    delx = (*inter_it).x[0] - itl_up->second->x[0];
                    dely = (*inter_it).x[1] - itl_up->second->x[1];
                    delz = (*inter_it).x[2] - itl_up->second->x[2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                fpair = pot->toForce(
                        atom_type::getTypeIdByType((*inter_it).type),
                                atom_type::getTypeIdByType(itl_up->second->type),
                                dist2, (*inter_it).df + itl_up->second->df);

                (*inter_it).f[0] += delx * fpair;
                (*inter_it).f[1] += dely * fpair;
                (*inter_it).f[2] += delz * fpair;

                        itl_up->second->f[0] -= delx * fpair;
                        itl_up->second->f[1] -= dely * fpair;
                        itl_up->second->f[2] -= delz * fpair;
                    }
                }
            }
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
    int xstart = p_domain->dbx_lattice_size_ghost[0];
    int ystart = p_domain->dbx_lattice_size_ghost[1];
    int zstart = p_domain->dbx_lattice_size_ghost[2];
    std::cout << "print_force" << std::endl;
    for (int k = zstart; k < p_domain->dbx_lattice_size_sub_box[2] + zstart; k++) {
        for (int j = ystart; j < p_domain->dbx_lattice_size_sub_box[1] + ystart; j++) {
            for (int i = xstart; i < p_domain->dbx_lattice_size_sub_box[0] + xstart; i++) {
                kk = atom_list->IndexOf3DIndex(i, j, k);
                AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                outfile << atom_.f[0] << " " << atom_.f[1] << " " << atom_.f[2] << std::endl;
            }
        }
    }
    outfile.close();
}

void atom::setv(int lat[4], double direction[3], double energy) {
    long kk;
    if ((lat[0] * 2) >= p_domain->dbx_lattice_coord_sub_box_region.x_low &&
        (lat[0] * 2) < (p_domain->dbx_lattice_coord_sub_box_region.x_low + p_domain->dbx_lattice_size_sub_box[0])
        && lat[1] >= p_domain->dbx_lattice_coord_sub_box_region.y_low &&
        lat[1] < (p_domain->dbx_lattice_coord_sub_box_region.y_low + p_domain->dbx_lattice_size_sub_box[1])
        && lat[2] >= p_domain->dbx_lattice_coord_sub_box_region.z_low &&
        lat[2] < (p_domain->dbx_lattice_coord_sub_box_region.z_low + p_domain->dbx_lattice_size_sub_box[2])) {
        kk = (atom_list->IndexOf3DIndex(lat[0] * 2 - p_domain->dbx_lattice_coord_ghost_region.x_low,
                                        lat[1] - p_domain->dbx_lattice_coord_ghost_region.y_low,
                                        lat[2] - p_domain->dbx_lattice_coord_ghost_region.z_low) + lat[3]);
        // todo verify the position.
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
        double v_ = sqrt(energy / atom_type::getAtomMass(atom_.type) / mvv2e); // the unit of v is 100m/s
        double d_ = sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
        atom_.v[0] += v_ * direction[0] / sqrt(d_);
        atom_.v[1] += v_ * direction[1] / sqrt(d_);
        atom_.v[2] += v_ * direction[2] / sqrt(d_);
    }
}

void atom::sendForce() {
    ForcePacker force_packer(getAtomListRef(), atom_list->sendlist, atom_list->recvlist);
    comm::neiSendReceive<double, MPI_DOUBLE>(&force_packer, MPIDomain::toCommProcess(), p_domain->rank_id_neighbours);
}
