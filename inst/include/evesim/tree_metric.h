//  SPDX-FileCopyrightText: 2024 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
//  SPDX-License-Identifier: MIT


#ifndef TREE_METRIC_H_INCLUDED
#define TREE_METRIC_H_INCLUDED
#pragma once

#include <type_traits>
#include "evesim/rutils.h"
#include "evesim/tree.h"


namespace tres_sim {

  namespace tree_metric {


    struct cophenetic {
      static void apply(const tree_t& tree, RcppParallel::RMatrix<double> D);
      Rcpp::NumericMatrix operator()(const tree_t& tree) const;

    private:
      struct dij_t {
        explicit dij_t(const tree_t& tree);

        // returns distance between tips i and j
        double operator()(double const* const Di, int i, int j) const;

        const tree_t& tree;
        std::vector<int> visitors;
      };
    };

    struct ed {
      static void apply(const tree_t& tree, RcppParallel::RVector<double> D);
      Rcpp::NumericVector operator()(const tree_t& tree) const;
    };


    struct nnd {
      static void apply(const tree_t& tree, RcppParallel::RVector<double> D);
      Rcpp::NumericVector operator()(const tree_t& tree) const;
    };


    struct pd {
      static void apply(const tree_t& tree, double& D);
      double operator()(const tree_t& tree) const;
    };


    struct mpd {
      static void apply(const tree_t& tree, double& D);
      double operator()(const tree_t& tree) const;
    };


    void set_dim_names(const Rcpp::RObject& MV, const tree_t& tree);
    inline void set_dim_names(const double, const tree_t&) {}


    template <typename Metric>
    inline auto metric(Metric&& m, const tree_t& tree) 
    -> std::invoke_result_t<Metric, const tree_t&> {
      auto ret = m(tree);
      set_dim_names(ret, tree);
      return ret;
    }


  } // namespace tree_metric

} // namespace tres_sim

#endif
