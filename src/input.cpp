#include <iostream>

#include <mpi.h>
#include <logs/logs.h>

#include "input.h"
#include "atom.h"

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
    readAtoms(fs, head, _atom, p_domain);
    fs.close();
}

inp_header input::readHeader(std::fstream &fs) {
    inp_header head;
    fs.read(reinterpret_cast<char *>(&head), sizeof(head));
    return head;
}

struct AtomForRead {
    long id;
    long type;
    double x[DIMENSION];
    double v[DIMENSION];
    double f[DIMENSION];
};

void input::readAtoms(std::fstream &fs, const inp_header head, atom *_atom, comm::BccDomain *p_domain) {
    // set all lattice site as invalid.
    _atom->atom_list->foreachSubBoxAtom([_atom](_type_atom_index gid) {
        MD_LOAD_ATOM_VAR(_atom_ref, _atom->atom_list, gid);
        MD_SET_ATOM_TYPE(_atom_ref, gid, atom_type::INVALID);
    });

    AtomForRead atom_bin;
    AtomElement atom;
    for (long i = 0; i < head.atom_count; i++) {
        fs.read(reinterpret_cast<char *>(&atom_bin), sizeof(AtomForRead));
        atom.id = atom.id;
        atom.type = atom.type;
        atom.x[0] = atom.x[0];
        atom.x[1] = atom.x[1];
        atom.x[2] = atom.x[2];
        atom.v[0] = atom.v[0];
        atom.v[1] = atom.v[1];
        atom.v[2] = atom.v[2];
        _atom->addAtom(p_domain, atom);
    }
}
