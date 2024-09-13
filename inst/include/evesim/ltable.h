//  SPDX-FileCopyrightText: 2024 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
//  SPDX-License-Identifier: MIT

#ifndef LTABLE_H_INCLUDED
#define LTABLE_H_INCLUDED
#pragma once

#include "Rcpp.h"
// [[Rcpp::plugins("cpp17")]]
#include "RcppParallel.h"
#include <vector>
#include <algorithm>


namespace tres_sim {
  

  template <bool descending>
  class ltable_view;

  
  // descending copy of R-ltable_t
  class ltable_t {
  public:
    struct entry_t {
      double t;
      int ancestor;
      double death;
      int label;      // signed self-index, one-based
    };

    ltable_t(const ltable_t&) = default;
    ltable_t(ltable_t&&) = default;
    ltable_t& operator=(const ltable_t&) = default;
    ltable_t& operator=(ltable_t&&) = default;

    explicit ltable_t(double age) : entries_{{age, -1, 0.0, 1}, {age, -1, 0.0, -2}} {}

    template <typename LTABLE_VIEW>
    ltable_t(const LTABLE_VIEW& LV) : entries_(LV.size()) {
      for (int i = 0; i < LV.size(); ++i) {
        entries_[i].t = LV.t(i);
        entries_[i].ancestor = LV.ancestor(i);
        entries_[i].death = LV.death(i);
        entries_[i].label = LV.label(i);
      }
      entries_[1].ancestor = 0;   // fix legacy quirk
    }

    int size() const noexcept { return static_cast<int>(entries_.size()); }
    double age() const noexcept { return entries_[0].t; }

    entry_t& operator[](int i) { assert((i >= 0) && (i < size())); return entries_[i]; }
    const entry_t& operator[](int i) const { assert((i >= 0) && (i < size())); return entries_[i]; }

  private:
    friend class sim_table_t;
    std::vector<entry_t> entries_;
  };


  // descending view into R-Ltable, zero-copy 
  template <bool descending>
  class ltable_view {
  public:
    static constexpr bool descending_source = descending;

    // precondition: conformant R ltable_t in LR
    ltable_view(const Rcpp::NumericMatrix& LR, double age);

    int size() const noexcept { return size_; }
    double age() const;

    // scalar accessors
    double t(int idx) const;
    int ancestor(int idx) const { return std::abs(static_cast<int>(LR_.column(1)[idx])) - 1; };
    double death(int idx) const;
    int label(int idx) const  { return (idx == 1) ? 2 : LR_.column(2)[idx]; };

    ltable_t::entry_t operator()(int idx) const {
      return { t(idx), ancestor(idx), death(idx), label(idx) };
    }

    // vector accessors to raw data
    auto raw_t() const noexcept { return LR_.column(0); }
    auto raw_ancestor() const noexcept { return LR_.column(1); }
    auto raw_death() const noexcept { return LR_.column(3); }
    auto raw_label() const noexcept { return LR_.column(2); }
        
    ltable_t to_ltable() const {
      return ltable_t(*this);
    }

    // returns descending L-table
    Rcpp::NumericMatrix to_matrix() const {
      auto LR = Rcpp::NumericMatrix(size(), 4);
      for (int i = 0; i < size(); ++i) {
        LR(i, 0) = t(i);
        LR(i, 1) = (LR_(i, 2) < 0) ? -LR_(i, 1) : +LR_(i, 1);   // signed ancestor
        LR(i, 2) = LR_(i, 2);   // signed self-index
        LR(i, 3) = (death(i) == 0.0) ? -1.0 : death(i);
      }
      if (LR_(0,0) != age()) {
        // was ascending - reverse times
        for (int i = 0; i < size(); ++i) {
          LR(i, 0) = age() - LR(i, 0);
          LR(i, 3) = (LR(i, 3) == -1.0) ? -1.0 : age() - LR(i, 3);
        }
      }
      return LR;
    }

  private:
    const RcppParallel::RMatrix<double> LR_;
    const double ofs_;    // age offset
    int size_;
  };


  // specializations for descending source
  template <>
  inline ltable_view<true>::ltable_view(const Rcpp::NumericMatrix& LR, double age) : LR_(LR), ofs_(age - LR(0,0)) {
    auto ct = LR_.column(0);
    auto it = std::lower_bound(ct.begin(), ct.end(), 0.0, [ofs = ofs_](double t, double val) {
      return (t + ofs) > val;
    });
    size_ = std::distance(ct.begin(), it);
  }

  template <> inline double ltable_view<true>::age() const { return LR_(0,0) + ofs_; };
  template <> inline double ltable_view<true>::t(int idx) const { return LR_(idx, 0) + ofs_; };
  template <> inline double ltable_view<true>::death(int idx) const { return std::max(0.0, LR_(idx, 3) + ofs_); }

  
  // specializations for ascending source
  template <>
  inline ltable_view<false>::ltable_view(const Rcpp::NumericMatrix& LR, double age) : LR_(LR), ofs_(age) {
    auto ct = LR_.column(0);
    auto it = std::lower_bound(ct.begin(), ct.end(), 0.0, [ofs = ofs_](double t, double val) {
      return (t - ofs) < val;
    });
    size_ = std::distance(ct.begin(), it);
  }

  template <> inline double ltable_view<false>::age() const { return ofs_; };
  template <> inline double ltable_view<false>::t(int idx) const { return ofs_ - LR_(idx, 0); };
  template <> inline double ltable_view<false>::death(int idx) const { 
    const auto d = LR_(idx, 3); 
    return ((d < 0.0) || (d > ofs_)) ? 0.0 : ofs_ - d;
  }


} // namespace tres_sim

#endif  // LTABLE_H_INCLUDED
