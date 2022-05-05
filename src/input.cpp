#include <iostream>

#include <mpi.h>
#include <logs/logs.h>

#include "input.h"
#include "atom.h"
#include "utils/mpi_data_types.h"

struct AtomForRead {
    _type_atom_id id;
    int64_t type; // padding int32 to int64.
    double x[DIMENSION];
    double v[DIMENSION];
    double f[DIMENSION];
};

input::input() {}

input::~input() {}

/**
 * This function will read atoms data from file into @param _atom.
 * The implementation does not guarantee the relationship of data in file and the simulation config (e.g. box size).
 * For example, if we read atoms data from a 40*40*40 box, but the simulation box is 80*80*80,
 * it may involve a undefined behavior.
 *
 * @param read_inp_path input file path containing all atoms.
 * @param _atom atoms container
 * @param p_domain the simulation domain
 */
void input::readPhaseSpace(const std::string read_inp_path, atom *_atom, comm::BccDomain *p_domain) {
    kiwi::logs::e("read-phase-space", "Opening phase space file {}.\n", read_inp_path);
    std::fstream fs(read_inp_path, std::ios::in | std::ios::binary);
    if (!fs.is_open()) {
        kiwi::logs::e("read-phase-space", "Could not open phaseSpaceFile {}.\n ", read_inp_path);
        MPI_Abort(MPI_COMM_WORLD, 1);
        return;
    }
    kiwi::logs::e("read-phase-space", "Reading phase space file {}.\n", read_inp_path);

    const inp_header head = readHeader(fs);
    assert(head.atom_size == sizeof(AtomForRead));
    const _type_atom_count atom_count = readAtoms(fs, head, _atom, p_domain);
    fs.close();

    // check results
    checkAtomRead(atom_count, p_domain);
}

void input::checkAtomRead(const _type_atom_count atom_added, comm::BccDomain *p_domain) {
    _type_atom_count total_atom_count = 0;
    MPI_Reduce(&atom_added, &total_atom_count, 1, mpi_types::mpi_type_atoms_count, MPI_SUM, MASTER_PROCESSOR,
               MPI_COMM_WORLD);

    const _type_atom_count expected_n_global_atoms =
            2 * p_domain->phase_space[0] * p_domain->phase_space[1] * p_domain->phase_space[2];
    if (kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR && total_atom_count != expected_n_global_atoms) {
        kiwi::logs::e(MASTER_PROCESSOR, "read-phase-space", "wrong atom number {} in the system, unexpect {}",
                      total_atom_count, expected_n_global_atoms);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

inp_header input::readHeader(std::fstream &fs) {
    inp_header head;
    fs.read(reinterpret_cast<char *>(&head), sizeof(head));
    return head;
}


_type_atom_count input::readAtoms(std::fstream &fs, const inp_header head, atom *_atom, comm::BccDomain *p_domain) {
    // set all lattice site as invalid.
    _atom->atom_list->foreachSubBoxAtom([_atom](_type_atom_index gid) {
        MD_LOAD_ATOM_VAR(_atom_ref, _atom->atom_list, gid);
        MD_SET_ATOM_TYPE(_atom_ref, gid, atom_type::INVALID);
    });

    AtomForRead atom_bin;
    AtomElement atom;
    _type_atom_count n = 0;
    for (long i = 0; i < head.atom_count; i++) {
        fs.read(reinterpret_cast<char *>(&atom_bin), sizeof(AtomForRead));
        atom.id = atom_bin.id;
        atom.type = static_cast<_type_atom_type_enum>(atom_bin.type);
        atom.x[0] = atom_bin.x[0];
        atom.x[1] = atom_bin.x[1];
        atom.x[2] = atom_bin.x[2];
        atom.v[0] = atom_bin.v[0];
        atom.v[1] = atom_bin.v[1];
        atom.v[2] = atom_bin.v[2];
        if (_atom->addAtom(p_domain, atom)) {
            n++;
        }
    }
    return n;
}
