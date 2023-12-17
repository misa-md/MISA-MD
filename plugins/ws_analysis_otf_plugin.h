//
// Created by genshen on 2023/12/12.
//

#ifndef MISA_MD_WS_ANALYSIS_OTF_PLUGIN_H
#define MISA_MD_WS_ANALYSIS_OTF_PLUGIN_H

#include <vector>

#include "atom/atom_set.h"
#include "lattice/ws_utils.h"
#include "plugin_api.h"

namespace plugins {
  static constexpr char TagNone = 0;
  static constexpr char TagCop = 1;

  /**
   * ws analysis on the fly.
   */
  class WSAnalysisOTFPlugin final : public IOPlugin {
  public:
    explicit WSAnalysisOTFPlugin(const BccLattice &_lat, comm::Domain *_p_domain) : lat(_lat), p_domain(_p_domain) {}

    ~WSAnalysisOTFPlugin() override = default;

    void init() override {
      const _type_atom_count n = lat._size_sub_box_x * lat._size_sub_box_y * lat._size_sub_box_z;
      set.resize(n);
      std::fill(set.begin(), set.end(), TagNone);
    }

    bool filter_atom(const _type_atom_location pos[3]) override {
      const auto near_atom_index = ws::getNearestLatIndexInSubBox(lat, pos, p_domain);
      if (near_atom_index == box::IndexNotExists) {
        return true;
      }
      if (set[near_atom_index] == TagNone) {
        set[near_atom_index] = TagCop;
        return false;
      } else {
        return true; // can output the atom
      }
      // todo set vacancy.
    }

  private:
    std::vector<char> set;
    const BccLattice &lat;
    comm::Domain *p_domain;
  };
} // namespace plugins

#endif // MISA_MD_WS_ANALYSIS_OTF_PLUGIN_H
