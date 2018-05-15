//
// Created by genshen on 5/7/18.
//

#include "atom_element.h"

bool AtomElement::isInterElement() const {
    return x[0] == COORDINATE_ATOM_OUT_BOX;
}
