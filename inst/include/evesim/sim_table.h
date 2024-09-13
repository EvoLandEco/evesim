//  SPDX-FileCopyrightText: 2024 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
//  SPDX-License-Identifier: MIT

#ifndef SIM_TABLE_H_INCLUDED
#define SIM_TABLE_H_INCLUDED
#pragma once

#include <array>
#include <vector>
#include "evesim/ltable.h"
#include "evesim/tree_metric.h"


namespace tres_sim {


  class sim_table_t {
  public:
    explicit sim_table_t(double age);
    explicit sim_table_t(const ltable_t& rhs);
    explicit sim_table_t(ltable_t&& rhs);

    // pune constructor
    sim_table_t(const sim_table_t& rhs, double age);
    ~sim_table_t() = default;

    // events
    void speciation(int ances, double t);
    void extinction(int specie, double t);

    double age() const noexcept { return ltable_[0].t; };
    auto nspecie() const noexcept { return specie_[0] + specie_[1]; };
    auto nclade_specie() const noexcept { return specie_; };
    int size() const noexcept { return ltable_.size(); }
  
    const tree_t& tree() const;
    tree_t full_tree() const;

    const ltable_t& ltable() const { return ltable_; };

    // one-based
    template <typename OIT> void specie_label(OIT oit) const {
      for (int i = 0; i < nspecie(); ++i) {
        *oit++ = tip_map_[i] + 1;
      }
    }

    const auto& tip_map() const noexcept { return tip_map_; }

  private:
    bool dirty() const noexcept { return cached_tree_.nodes.empty(); }

    ltable_t ltable_;   // ascending ltable
    mutable tree_t cached_tree_;  // cached drop extinct tree
    std::vector<int> tip_map_;  // extant specie tip label (unsorted)
    std::array<int, 2> specie_ = { 0, 0 };  // specie species per root-clade
  };


}

#endif  // SIM_TABLE_H_INCLUDED
