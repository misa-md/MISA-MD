#include <iostream>
#include <iterator>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <utils/mpi_domain.h>
#include <logs/logs.h>
#include <eam.h>
#include <comm/comm.hpp>

#include "atom.h"
#include "lattice/ws_utils.h"
#include "pack/rho_packer.h"
#include "pack/force_packer.h"
#include "pack/df_embed_packer.h"
#include "arch/hardware_accelerate.hpp" // use hardware(eg.GPU, MIC,Sunway slave cores.) to achieve calculate accelerating.

atom::atom(comm::BccDomain *domain)
        : AtomSet(domain->lattice_const * domain->cutoff_radius_factor,
                  domain->ghost_extended_lattice_size,
                  domain->sub_box_lattice_size,
                  domain->lattice_size_ghost),
          p_domain(domain) {
    if (isAccelerateSupport()) {
        accelerateInit(p_domain->dbx_sub_box_lattice_region.x_low,
                       p_domain->dbx_sub_box_lattice_region.y_low,
                       p_domain->dbx_sub_box_lattice_region.z_low,
                       p_domain->dbx_sub_box_lattice_size[0],
                       p_domain->dbx_sub_box_lattice_size[1],
                       p_domain->dbx_sub_box_lattice_size[2],
                       p_domain->dbx_ghost_ext_lattice_region.x_low,
                       p_domain->dbx_ghost_ext_lattice_region.y_low,
                       p_domain->dbx_ghost_ext_lattice_region.z_low,
                       p_domain->dbx_ghost_extended_lattice_size[0],
                       p_domain->dbx_ghost_extended_lattice_size[1],
                       p_domain->dbx_ghost_extended_lattice_size[2]);
    }
}

int atom::decide() {
    inter_atom_list->clearGhost();
    int nflag = 0;
    double dist;
    double xtemp, ytemp, ztemp;

    //对本地晶格点原子进行判断，看是否运动为间隙原子
    for (_type_atom_index k = 0; k < p_domain->dbx_sub_box_lattice_size[2]; k++) {
        for (_type_atom_index j = 0; j < p_domain->dbx_sub_box_lattice_size[1]; j++) {
            for (_type_atom_index i = 0; i < p_domain->dbx_sub_box_lattice_size[0]; i++) {
//                kk = atom_list->IndexOf3DIndex(i, j, k);
                AtomElement &atom_ = atom_list->getAtomEleBySubBoxIndex(i, j, k); // todo long type
                if (!atom_.isInterElement()) {
                    xtemp = (i + p_domain->dbx_sub_box_lattice_region.x_low) * 0.5 * p_domain->lattice_const;
                    ytemp = (j + p_domain->dbx_sub_box_lattice_region.y_low + (i % 2) * 0.5) *
                            p_domain->lattice_const;
                    ztemp = (k + p_domain->dbx_sub_box_lattice_region.z_low + (i % 2) * 0.5) *
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
    for (_type_atom_index i = 0; i < atom_list->cap(); i++) {
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

void atom::computeEam(eam *pot, double &comm) {
    double starttime, stoptime;
    inter_atom_list->makeIndex(atom_list, p_domain); // create index for inter atom and inter ghost atoms.

    latRho(pot);
    interRho(pot);
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
        comm::neiSendReceive<double>(&rho_packer, MPIDomain::toCommProcess(),
                                     MPI_DOUBLE, p_domain->rank_id_neighbours, true);
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
    latDf(pot);

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
        comm::neiSendReceive<double>(&packer, MPIDomain::toCommProcess(), MPI_DOUBLE, p_domain->rank_id_neighbours);
        stoptime = MPI_Wtime();
        comm += stoptime - starttime;
    }

    // force for local lattice.
    latForce(pot);

    /*sprintf(tmp, "f.atom");  // 2.todo remove start.
      outfile.open(tmp);
      for(int i = 0; i < phi_spline->n; i++){
         outfile << i << " " << phi_spline->spline[i][6] << std::endl;
      }
      outfile.close();*/ // 2.todo remove end.
    //间隙原子计算嵌入能和对势带来的力
    interForce(pot);

    // send force
    starttime = MPI_Wtime();
    ForcePacker force_packer(getAtomListRef(), atom_list->sendlist, atom_list->recvlist);
    comm::neiSendReceive<double>(&force_packer, MPIDomain::toCommProcess(),
                                 MPI_DOUBLE, p_domain->rank_id_neighbours, true);
    stoptime = MPI_Wtime();
    comm += stoptime - starttime;
}

void atom::latRho(eam *pot) {
    double xtemp, ytemp, ztemp;
    double delx, dely, delz;
    _type_atom_index kk;
    double dist2;
    int xstart = p_domain->dbx_lattice_size_ghost[0];
    int ystart = p_domain->dbx_lattice_size_ghost[1];
    int zstart = p_domain->dbx_lattice_size_ghost[2];

    // 本地晶格点上的原子计算电子云密度
    if (isAccelerateSupport()) {
//     fixme  accelerateEamRhoCalc(&(rho_spline->n), atom_list, &_cutoffRadius,
//                             &(rho_spline->invDx), rho_spline->values); // fixme
    } else { // calculate electron density use cpu only.
        for (int k = zstart; k < p_domain->dbx_sub_box_lattice_size[2] + zstart; k++) {
            for (int j = ystart; j < p_domain->dbx_sub_box_lattice_size[1] + ystart; j++) {
                for (int i = xstart; i < p_domain->dbx_sub_box_lattice_size[0] + xstart; i++) {
                    kk = atom_list->lattice.IndexOf3DIndex(i, j, k);
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
                                atom_central.rho += pot->chargeDensity(
                                        atom_type::getTypeIdByType(atom_neighbour.type), dist2);
                                atom_neighbour.rho += pot->chargeDensity(
                                        atom_type::getTypeIdByType(atom_central.type), dist2);
                            }
                        }
                    }
                }
            }
        }
    }
}

void atom::interRho(eam *pot) {
    double delx, dely, delz;
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
            inter_it->rho += pot->chargeDensity(atom_type::getTypeIdByType(atom_near.type), dist2);
            atom_near.rho += pot->chargeDensity(atom_type::getTypeIdByType((*inter_it).type), dist2);
        }

        _type_atom_index x, y, z;
        atom_list->lattice.get3DIndexByLinearIndex(near_atom_index, x, y, z);
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
                    (*inter_it).rho += pot->chargeDensity(
                            atom_type::getTypeIdByType(lat_nei_atom.type), dist2);
                    lat_nei_atom.rho += pot->chargeDensity(
                            atom_type::getTypeIdByType((*inter_it).type), dist2);
                }
            }
        }

        // rho contribution of neighbour atoms in the same bucket.
        {
            const inter_map_range inter_map_range = inter_atom_list->inter_map.equal_range(near_atom_index);
            for (inter_map_range_itl bucket_nei_itl = inter_map_range.first;
                 bucket_nei_itl != inter_map_range.second; ++bucket_nei_itl) {
                if (bucket_nei_itl->second->id == inter_it->id) {
                    continue; // can not be itself.
                }
                delx = (*inter_it).x[0] - bucket_nei_itl->second->x[0];
                dely = (*inter_it).x[1] - bucket_nei_itl->second->x[1];
                delz = (*inter_it).x[2] - bucket_nei_itl->second->x[2];
                dist2 = delx * delx + dely * dely + delz * delz;
                if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                    (*inter_it).rho += pot->chargeDensity(
                            atom_type::getTypeIdByType(bucket_nei_itl->second->type), dist2);
                }
            }
        }

        // rho between inter atoms and inter atoms (use full neighbour index).
        AtomNei::iterator nei_half_itl_end = neighbours->end(false, x, y, z);
        for (AtomNei::iterator nei_itl = neighbours->begin(false, x, y, z);
             nei_itl != nei_half_itl_end; ++nei_itl) {
            const _type_atom_index inter_nei_id = atom_list->lattice.IndexOf3DIndex(
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
                    (*inter_it).rho += pot->chargeDensity(atom_type::getTypeIdByType(itl->second->type), dist2);
//                    itl->second->rho += pot->rhoContribution(atom_type::getTypeIdByType((*inter_it).type), dist2);
                }
            }
        }
        //计算间隙原子嵌入能导数
        dfEmbed = pot->dEmbedEnergy(atom_type::getTypeIdByType((*inter_it).type), (*inter_it).rho);
        (*inter_it).df = dfEmbed;
    }
}

void atom::latDf(eam *pot) {
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
        for (int k = zstart; k < p_domain->dbx_sub_box_lattice_size[2] + zstart; k++) {
            for (int j = ystart; j < p_domain->dbx_sub_box_lattice_size[1] + ystart; j++) {
                for (int i = xstart; i < p_domain->dbx_sub_box_lattice_size[0] + xstart; i++) {
                    kk = atom_list->lattice.IndexOf3DIndex(i, j, k);
                    AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                    if (atom_.isInterElement()) {
                        continue;
                    }
                    dfEmbed = pot->dEmbedEnergy(atom_type::getTypeIdByType(atom_.type), atom_.rho); // fixme
                    atom_.df = dfEmbed;
                }
            }
        }
    }
}

void atom::latForce(eam *pot) {
    double xtemp, ytemp, ztemp;
    double delx, dely, delz;
    _type_atom_index kk;
    double dist2;
    double fpair;

    int xstart = p_domain->dbx_lattice_size_ghost[0];
    int ystart = p_domain->dbx_lattice_size_ghost[1];
    int zstart = p_domain->dbx_lattice_size_ghost[2];

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

        for (int k = zstart; k < p_domain->dbx_sub_box_lattice_size[2] + zstart; k++) {
            for (int j = ystart; j < p_domain->dbx_sub_box_lattice_size[1] + ystart; j++) {
                for (int i = xstart; i < p_domain->dbx_sub_box_lattice_size[0] + xstart; i++) {
                    kk = atom_list->lattice.IndexOf3DIndex(i, j, k);
                    AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
                    xtemp = atom_.x[0];
                    ytemp = atom_.x[1];
                    ztemp = atom_.x[2];
                    if (atom_.isInterElement()) {
                        continue;
                    }

                    // force between lattice atoms and lattice atoms.
                    AtomNei::iterator nei_itl_end = neighbours->end(true, i, j, k);
                    for (AtomNei::iterator nei_itl = neighbours->begin(true, i, j, k);
                         nei_itl != nei_itl_end; ++nei_itl) {
                        AtomElement &atom_n = *nei_itl;
                        delx = xtemp - atom_n.x[0];
                        dely = ytemp - atom_n.x[1];
                        delz = ztemp - atom_n.x[2];
                        dist2 = delx * delx + dely * dely + delz * delz;
                        if (dist2 < (_cutoffRadius * _cutoffRadius) && !atom_n.isInterElement()) {
                            fpair = pot->toForce(atom_type::getTypeIdByType(atom_.type),
                                                 atom_type::getTypeIdByType(atom_n.type),
                                                 dist2, atom_.df, atom_n.df);

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

void atom::interForce(eam *pot) {
    double delx, dely, delz;
    double dist2;
    double fpair;

    //间隙原子计算嵌入能和对势带来的力
    _type_atom_index _atom_near_index;
    for (_type_inter_list::iterator inter_it = inter_atom_list->inter_list.begin();
         inter_it != inter_atom_list->inter_list.end(); inter_it++) {
        _atom_near_index = ws::findNearLatIndexInSubBox(atom_list, *inter_it, p_domain);
        if (_atom_near_index == box::IndexNotExists) {
            assert(false);
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
            fpair = pot->toForce(
                    atom_type::getTypeIdByType((*inter_it).type),
                    atom_type::getTypeIdByType(atom_central.type),
                    dist2, (*inter_it).df, atom_central.df);

            (*inter_it).f[0] += delx * fpair;
            (*inter_it).f[1] += dely * fpair;
            (*inter_it).f[2] += delz * fpair;

            atom_central.f[0] -= delx * fpair;
            atom_central.f[1] -= dely * fpair;
            atom_central.f[2] -= delz * fpair;
        }

        _type_atom_index x, y, z;
        atom_list->lattice.get3DIndexByLinearIndex(_atom_near_index, x, y, z);
        // force between inter atoms and lattice atoms (use full neighbour index).
        AtomNei::iterator nei_full_itl_end = neighbours->end(false, x, y, z);
        for (AtomNei::iterator nei_itl = neighbours->begin(false, x, y, z);
             nei_itl != nei_full_itl_end; ++nei_itl) {
            AtomElement &lattice_neighbour = *nei_itl; // this is a lattice atom.
            delx = (*inter_it).x[0] - lattice_neighbour.x[0];
            dely = (*inter_it).x[1] - lattice_neighbour.x[1];
            delz = (*inter_it).x[2] - lattice_neighbour.x[2];
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius) && !lattice_neighbour.isInterElement()) {
                fpair = pot->toForce(
                        atom_type::getTypeIdByType((*inter_it).type),
                        atom_type::getTypeIdByType(lattice_neighbour.type),
                        dist2, (*inter_it).df, lattice_neighbour.df);

                (*inter_it).f[0] += delx * fpair;
                (*inter_it).f[1] += dely * fpair;
                (*inter_it).f[2] += delz * fpair;

                lattice_neighbour.f[0] -= delx * fpair;
                lattice_neighbour.f[1] -= dely * fpair;
                lattice_neighbour.f[2] -= delz * fpair;
            }
        }
        // force contribution of neighbour atoms in the same bucket.
        {
            const inter_map_range inter_map_range = inter_atom_list->inter_map.equal_range(_atom_near_index);
            for (inter_map_range_itl bucket_nei_itl = inter_map_range.first;
                 bucket_nei_itl != inter_map_range.second; ++bucket_nei_itl) {
                if (bucket_nei_itl->second->id == inter_it->id) {
                    continue; // can not be itself.
                }
                delx = (*inter_it).x[0] - bucket_nei_itl->second->x[0];
                dely = (*inter_it).x[1] - bucket_nei_itl->second->x[1];
                delz = (*inter_it).x[2] - bucket_nei_itl->second->x[2];
                dist2 = delx * delx + dely * dely + delz * delz;
                if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                    fpair = pot->toForce(
                            atom_type::getTypeIdByType((*inter_it).type),
                            atom_type::getTypeIdByType(bucket_nei_itl->second->type),
                            dist2, (*inter_it).df, bucket_nei_itl->second->df);

                    (*inter_it).f[0] += delx * fpair;
                    (*inter_it).f[1] += dely * fpair;
                    (*inter_it).f[2] += delz * fpair;
                }
            }
        }

        // force between inter atoms and inter atoms(including inter ghost atom) (use full neighbour index).
        AtomNei::iterator nei_half_itl_end = neighbours->end(false, x, y, z);
        for (AtomNei::iterator nei_itl = neighbours->begin(false, x, y, z);
             nei_itl != nei_half_itl_end; ++nei_itl) {
            const _type_atom_index inter_nei_id = atom_list->lattice.IndexOf3DIndex(
                    nei_itl.cur_index_x, nei_itl.cur_index_y,
                    nei_itl.cur_index_z); // get index of the neighbour lattice.
            inter_map_range inter_map_range_up = inter_atom_list->inter_map.equal_range(inter_nei_id);
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
                            dist2, (*inter_it).df, itl_up->second->df);

                    (*inter_it).f[0] += delx * fpair;
                    (*inter_it).f[1] += dely * fpair;
                    (*inter_it).f[2] += delz * fpair;
                }
            }
        }
    }
}

void atom::setv(const _type_lattice_coord lat[4], const double direction[3], const double energy) {
    long kk;
    if ((lat[0] * 2) >= p_domain->dbx_sub_box_lattice_region.x_low &&
        (lat[0] * 2) < (p_domain->dbx_sub_box_lattice_region.x_low + p_domain->dbx_sub_box_lattice_size[0])
        && lat[1] >= p_domain->dbx_sub_box_lattice_region.y_low &&
        lat[1] < (p_domain->dbx_sub_box_lattice_region.y_low + p_domain->dbx_sub_box_lattice_size[1])
        && lat[2] >= p_domain->dbx_sub_box_lattice_region.z_low &&
        lat[2] < (p_domain->dbx_sub_box_lattice_region.z_low + p_domain->dbx_sub_box_lattice_size[2])) {
        kk = (atom_list->lattice.IndexOf3DIndex(lat[0] * 2 - p_domain->dbx_ghost_ext_lattice_region.x_low,
                                                lat[1] - p_domain->dbx_ghost_ext_lattice_region.y_low,
                                                lat[2] - p_domain->dbx_ghost_ext_lattice_region.z_low) + lat[3]);
        // todo verify the position.
        AtomElement &atom_ = atom_list->getAtomEleByLinearIndex(kk);
        const double v_ = sqrt(2 * energy / atom_type::getAtomMass(atom_.type) / mvv2e); // the unit of v is A/ps (or 100m/s)
        const double d_ = sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
        atom_.v[0] += v_ * direction[0] / d_;
        atom_.v[1] += v_ * direction[1] / d_;
        atom_.v[2] += v_ * direction[2] / d_;
    }
}
