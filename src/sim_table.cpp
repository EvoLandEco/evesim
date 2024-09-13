//  SPDX-FileCopyrightText: 2024 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
//  SPDX-License-Identifier: MIT

#include <tuple>
#include "evesim/sim_table.h"


namespace tres_sim {


  sim_table_t::sim_table_t(double age) : ltable_(age), tip_map_{0,1}, specie_{1,1} {};


  sim_table_t::sim_table_t(const ltable_t& rhs) : sim_table_t(ltable_t(rhs)) {}


  sim_table_t::sim_table_t(ltable_t&& rhs) : ltable_(std::move(rhs)) {
    for (int i = 0; i < ltable_.size(); ++i) {
      const auto clade_specie = static_cast<int>(ltable_[i].label < 0);
      ++specie_[clade_specie];
      if (ltable_.entries_[i].death == 0.0) {
        tip_map_.push_back(i);
      }
    }
  }


  // prunning ctor
  sim_table_t::sim_table_t(const sim_table_t& rhs, double age) : ltable_(0.0) {
    ltable_.entries_.clear();
    const auto& ltable = rhs.ltable();
    for (int i = 0; i < ltable.size(); ++i) {
      const auto e = ltable.entries_[i];
      if (e.t > age) break;
      auto& ne = ltable_.entries_.emplace_back(e);
      if ((ne.death == 0.0) || (ne.death >= age)) {
        ne.death = 0.0;
        tip_map_.push_back(i);
        ++specie_[static_cast<int>(ne.label < 0)];
      }
    }
  }


  void sim_table_t::speciation(int specie, double t) {
    const auto ances = tip_map_[specie];
    tip_map_.push_back(ltable_.size());
    const auto clade = static_cast<int>(ltable_[ances].label < 0);
    const auto label = clade ? -(ltable_.size() + 1) : (ltable_.size() + 1);
    ltable_.entries_.push_back({ age() - t, ances, 0.0, label });
    ++specie_[clade];
    cached_tree_.clear();
  }


  void sim_table_t::extinction(int specie, double t) {
    auto pivot = tip_map_[specie];
    ltable_[pivot].death = age() - t;
    tip_map_.erase(tip_map_.begin() + specie);
    const auto clade = static_cast<int>(ltable_[pivot].label < 0);
    --specie_[clade];
    cached_tree_.clear();
  }


  const tree_t& sim_table_t::tree() const {
    if (dirty()) {
      cached_tree_ = tree_t::from(ltable_, true);
    }
    return cached_tree_;
  }


  tree_t sim_table_t::full_tree() const {
    return tree_t::from(ltable_, false);
  }


} // namespace tres_sim
