//
// Created by genshen on 2019-07-25.
//

#include "buffered_io.h"

BufferedFileWriter::BufferedFileWriter(kiwi::LocalStorage *p_local_storage, const unsigned long buffer_size)
        : p_local_storage(p_local_storage), buffer_size(buffer_size),
          next_buffer_index(0), total_atoms_written(0) {
    buffer = new atom_dump::AtomInfoDump[buffer_size];
    // todo get buffer size from local_storage's block size.
}

BufferedFileWriter::~BufferedFileWriter() {
    delete buffer;
}

void BufferedFileWriter::write(AtomElement *atom, const unsigned long time_step) {
    buffer[next_buffer_index].id = atom->id;
    buffer[next_buffer_index].step = time_step;
    buffer[next_buffer_index].type = atom->type;
    buffer[next_buffer_index].inter_type = 0; // normal
    buffer[next_buffer_index].atom_location[0] = atom->x[0];
    buffer[next_buffer_index].atom_location[1] = atom->x[1];
    buffer[next_buffer_index].atom_location[2] = atom->x[2];
    buffer[next_buffer_index].atom_velocity[0] = atom->v[0];
    buffer[next_buffer_index].atom_velocity[1] = atom->v[1];
    buffer[next_buffer_index].atom_velocity[2] = atom->v[2];
    next_buffer_index++;
    if (next_buffer_index == buffer_size) {
        // write data and reset index
        total_atoms_written += buffer_size;
        p_local_storage->writer.write(buffer, buffer_size * sizeof(atom_dump::AtomInfoDump));
        next_buffer_index = 0;
    }
}

void BufferedFileWriter::flush() {
    if (next_buffer_index == 0) {
        return;
    }
    total_atoms_written += next_buffer_index;
    p_local_storage->writer.write(buffer, next_buffer_index * sizeof(atom_dump::AtomInfoDump));
    next_buffer_index = 0;
}
