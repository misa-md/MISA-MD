//
// Created by genshen on 5/31/18.
//

#ifndef MISA_MD_ATOM_DUMP_TYPES_H
#define MISA_MD_ATOM_DUMP_TYPES_H


#include <cstddef>
#include <array>
#include <mpi.h>

#include "types/pre_define.h"

namespace atom_dump {
    typedef unsigned int type_dump_mask;
    // the dump_mask of dump option, it is associated with MPI_DataTypes array (as array index)
    constexpr type_dump_mask WithPositionMask = 1 << 0;
    constexpr type_dump_mask WithVelocityMask = 1 << 1;
    constexpr type_dump_mask WithForceMask = 1 << 2;;

    struct GlobalMetaData {
        size_t self_size; // in bytes
        size_t frame_meta_size; // in bytes
        size_t block_atoms; // atoms number in one block
        size_t atoms_num; // total atoms number
        size_t atom_item_bytes; // bytes of one atom element
        size_t mpi_ranks; // MPI ranks for running the simulation
        type_dump_mask mask; // dump mask (with force/velocity/position)
        unsigned int format_version;
        size_t global_header_size; // global header size in bytes.
        size_t local_size; // local header size in bytes.
        unsigned int frames; // frames in this file
    };

    struct FrameMetaData {
        size_t atoms_num; // atom count in current frame
        size_t atoms_num_hash_collision; // atoms number of hash collision
        size_t step; // step of current frame
        double time; // time of current frame
    };

    struct AtomInfoDump {
    public:
        _type_atom_id id;
        _type_atom_type type;
        _type_inter_type inter_type;
        _type_atom_location atom_location[DIMENSION]; // atom location
        _type_atom_velocity atom_velocity[DIMENSION]; // atom velocity
        _type_atom_force atom_force[DIMENSION]; // atom force
    };

    size_t getAtomDumpDataSize(const type_dump_mask mask);

    MPI_Datatype registerAtomDumpMPIDataType(const type_dump_mask mask);

    MPI_Datatype registerFrameMetaMPIDataType();

    static MPI_Datatype mpi_frame_meta_datatype = MPI_DATATYPE_NULL;
    static std::array<MPI_Datatype, 8> mpi_dump_types = {
            MPI_DATATYPE_NULL,
            MPI_DATATYPE_NULL,
            MPI_DATATYPE_NULL,
            MPI_DATATYPE_NULL,
            MPI_DATATYPE_NULL,
            MPI_DATATYPE_NULL,
            MPI_DATATYPE_NULL,
            MPI_DATATYPE_NULL,
    };
}


#endif //MISA_MD_ATOM_DUMP_TYPES_H
