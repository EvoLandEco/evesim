//  SPDX-FileCopyrightText: 2024 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
//  SPDX-License-Identifier: MIT

#include <vector>
#include "tres_sim/phylo.h"
#include "tres_sim/ltable.h"
#include "tres_sim/tree.h"
#include "tres_sim/tree_visit.h"


namespace tres_sim {

  namespace detail {

    Rcpp::List create_ape_phylo(size_t nnode) {
      auto phy = Rcpp::List::create(
        Rcpp::Named("edge") = Rcpp::IntegerMatrix(2 * nnode, 2),
        Rcpp::Named("edge.length") = Rcpp::NumericVector(2 * nnode),
        Rcpp::Named("Nnode") = static_cast<int>(nnode),
        Rcpp::Named("tip.label") = Rcpp::IntegerVector(nnode + 1),
        Rcpp::Named("root.edge") = 0
      );
      phy.attr("class") = "phylo";
      phy.attr("order") = "cladewise";
      return phy;
    }

  }


  phylo_t::phylo_t(Rcpp::List&& phy) :
    phy_(std::move(phy)),
    edge_(Rcpp::as<Rcpp::IntegerMatrix>(phy_["edge"])),
    edge_length_(Rcpp::as<Rcpp::NumericVector>(phy_["edge.length"])),
    tip_label_(Rcpp::as<Rcpp::IntegerVector>(phy_["tip.label"]))
  {}


  phylo_t::phylo_t(const tree_t& tree) : phylo_t(detail::create_ape_phylo(tree.nnode())) {
    auto& nodes = tree.nodes;
    if (nodes.empty()) return; 
    auto label_map = std::vector<int>(nodes.size(), -1);   // index map nodes -> label
    auto inode = tree.nnode() + 1;      // current inode label
    int ileaf = 0;                      // current leaf label
    int next_node = 0;
    auto e0 = this->e0();
    auto e1 = this->e1();
    auto len = this->len();
    auto tip_label = this->tip_label();
    visit(visit_order::clade_wise{}, tree, [&](int idx){
      auto& node = tree.nodes[idx];
      auto label = label_map[node.ances]; 
      if (label < 0) label = label_map[node.ances] = inode++;
      e0[next_node] = label + 1;   // one based
      e1[next_node] = idx;
      len[next_node] = nodes[node.ances].t - node.t;
      ++next_node;
      if (node.is_leaf()) {
        tip_label[ileaf] = std::abs(node.label);
        label_map[idx] = ileaf++;
      }
    });
    // remap e1 'labels'
    for (auto& e : e1) {
      e = label_map[e] + 1;   // one based
    }
  }


  Rcpp::List phylo_t::unwrap() {
    return std::move(phy_);
  }

}