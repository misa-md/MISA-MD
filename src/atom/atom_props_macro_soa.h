//
// Created by chu genshen on 2021/8/25.
//

#ifndef MISA_MD_ATOM_PROPS_MACRO_SOA_H
#define MISA_MD_ATOM_PROPS_MACRO_SOA_H

#include "../types/pre_define.h"
#include "atom_element.h"
#include "atom_prop_list.hpp"
#include "types/atom_types.h"

#define MD_HASH_LIST_DECLARE() \
AtomPropList<_type_atom_id> _atoms; /* atom id list */ \
AtomPropList<_type_atom_type_enum> atom_types; /* atom types list */ \
AtomPropList<_type_atom_location[3]> atom_x; /* atom position list */ \
AtomPropList<_type_atom_velocity[3]> atom_v; /* atom velocity list */ \
AtomPropList<_type_atom_force[3]> atom_f; /* atom force list */ \
AtomPropList<_type_atom_rho> atom_rho; /* atom rho list */ \
AtomPropList<_type_atom_df> atom_df; /* atom df list */

#define MD_HASH_LIST_INIT(lattice) \
_atoms(AtomPropList<_type_atom_id>(lattice)), \
atom_types(AtomPropList<_type_atom_type_enum>(lattice)), \
atom_x(AtomPropList<_type_atom_location[3]>(lattice)), \
atom_v(AtomPropList<_type_atom_velocity[3]>(lattice)), \
atom_f(AtomPropList<_type_atom_force[3]>(lattice)), \
atom_rho(AtomPropList<_type_atom_rho>(lattice)), \
atom_df(AtomPropList<_type_atom_df>(lattice))

#define MD_HASH_LIST_DESTROY(lattice) \
_atoms.destroyPropList(); \
atom_types.destroyPropList(); \
atom_x.destroyPropList(); \
atom_v.destroyPropList(); \
atom_f.destroyPropList(); \
atom_rho.destroyPropList(); \
atom_df.destroyPropList();

/**
 * id: global id
 * list: atom list pointer.
 */
#define MD_LOAD_ATOM_VAR(name, list, id) \
AtomPropList<_type_atom_id> &name ## _id = list->_atoms; \
AtomPropList<_type_atom_type_enum> &name ## _tp = list->atom_types; \
AtomPropList<_type_atom_location[3]> &name ## _x = list->atom_x; \
AtomPropList<_type_atom_velocity[3]> &name ## _v = list->atom_v; \
AtomPropList<_type_atom_force[3]> &name ## _f = list->atom_f; \
AtomPropList<_type_atom_rho> &name ## _rho = list->atom_rho; \
AtomPropList<_type_atom_df> &name ## _df = list->atom_df;

#define MD_GET_ATOM_ID(name, gid) \
name ## _id[gid]

#define MD_SET_ATOM_ID(name, gid, __id) \
name ## _id[gid] = __id

#define MD_GET_ATOM_TYPE(name, gid) \
name ## _tp[gid]

#define MD_SET_ATOM_TYPE(name, gid, __tp) \
name ## _tp[gid] = __tp

#define MD_IS_ATOM_TYPE_INTER(name, gid) \
(name ## _tp[gid] == atom_type::INVALID)

#define MD_GET_ATOM_X_ALL(name, gid) \
name ## _x[gid]

#define MD_GET_ATOM_X(name, gid, i) \
name ## _x[gid][i]

#define MD_SET_ATOM_X(name, gid, i, __x) \
name ## _x[gid][i] = __x

#define MD_ADD_ATOM_X(name, gid, i, __x) \
name ## _x[gid][i] += __x

#define MD_GET_ATOM_V(name, gid, i) \
name ## _v[gid][i]

#define MD_SET_ATOM_V(name, gid, i, __v) \
name ## _v[gid][i] = __v

#define MD_ADD_ATOM_V(name, gid, i, __v) \
name ## _v[gid][i] += __v

#define MD_GET_ATOM_F(name, gid, i) \
name ## _f[gid][i]

#define MD_SET_ATOM_F(name, gid, i, __f) \
name ## _f[gid][i] = __f

#define MD_ADD_ATOM_F(name, gid, i, __f) \
name ## _f[gid][i] += __f

#define MD_GET_ATOM_RHO(name, gid) \
name ## _rho[gid]

#define MD_SET_ATOM_RHO(name, gid, __rho) \
name ## _rho[gid] = __rho

#define MD_ADD_ATOM_RHO(name, gid, __rho) \
name ## _rho[gid] += __rho

#define MD_GET_ATOM_DF(name, gid) \
name ## _df[gid]

#define MD_SET_ATOM_DF(name, gid, __df) \
name ## _df[gid] = __df

inline AtomElement
convert_list_to_atom_element(const size_t gid, AtomPropList<_type_atom_id> &id_list,
                             AtomPropList<_type_atom_type_enum> &type_list,
                             AtomPropList<_type_atom_location[DIMENSION]> &location_list,
                             AtomPropList<_type_atom_velocity[DIMENSION]> &velocity_list,
                             AtomPropList<_type_atom_force[DIMENSION]> &force_list,
                             AtomPropList<_type_atom_rho> &rho_list, AtomPropList<_type_atom_df> &df_list) {
    AtomElement ele = AtomElement{.id = id_list[gid], .type = type_list[gid],
            .rho = rho_list[gid], .df = df_list[gid]};
    ele.x[0] = location_list[gid][0];
    ele.x[1] = location_list[gid][1];
    ele.x[2] = location_list[gid][2];

    ele.v[0] = velocity_list[gid][0];
    ele.v[1] = velocity_list[gid][1];
    ele.v[2] = velocity_list[gid][2];

    ele.f[0] = force_list[gid][0];
    ele.f[1] = force_list[gid][1];
    ele.f[2] = force_list[gid][2];

    return ele;
}

#define MD_TO_ATOM_ELEMENT(name, gid) \
convert_list_to_atom_element(gid, name ## _id, name ## _tp, name ## _x,  name ## _v, name ## _f, name ## _rho, name ## _df)

#endif //MISA_MD_ATOM_PROPS_MACRO_SOA_H
