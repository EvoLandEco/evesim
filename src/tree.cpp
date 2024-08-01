//  SPDX-FileCopyrightText: 2024 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
//  SPDX-License-Identifier: MIT

#include <utility>
#include "tres_sim/tree.h"
#include "tres_sim/ltable.h"
#include "tres_sim/phylo.h"


namespace tres_sim {

  namespace detail {

    tree_t to_tree(const ltable_t& ltable, std::true_type /*drop_extinct*/) {
      auto last_child = std::vector<int>(ltable.size(), 0);
      int tips = (ltable[0].death == 0.0) + (ltable[1].death == 0.0);
      for (int i = ltable.size() - 1; i > 1; --i) {
        const auto is_tip = ltable[i].death == 0.0;
        tips += is_tip;
        const auto ancestor = ltable[i].ancestor;
        if ((last_child[ancestor] == 0) && (is_tip || last_child[i])) {
          last_child[ancestor] = i;   // child with descendants or extant leaf
        }
      }
      auto nodes = std::vector<node_t>(2 * tips - 1);
      nodes[tips] = { ltable[0].t, -1, {0, 1}, 0 };    // a.k.a. root
      nodes[0] = { ltable[0].death, tips, {-1,-1}, ltable[0].label };
      nodes[1] = { ltable[1].death, tips, {-1,-1}, ltable[1].label };
      auto node_map = std::vector<int>(ltable.size());  // LL[LL[i].ances].t == nodes[node_map[i]].t
      node_map[0] = 0;
      node_map[1] = 1;
      int next_tip = 2;
      int next_node = tips + 1;
      for (int i = 2; i < ltable.size(); ++i) {
        if (const auto death = ltable[i].death; (death == 0.0) || last_child[i]) {
          // i has descendants or is extant leaf
          const auto ancestor = ltable[i].ancestor;
          auto const label = ltable[i].label;
          const auto ances = node_map[ancestor];
          if ((i == last_child[ancestor]) && (nodes[ances].t != 0.0)) {
            nodes[ances].t = death;
            nodes[ances].label = label;
            node_map[i] = ances;
          }
          else {
            const auto aa = nodes[ances].ances;
            nodes[next_node] = { ltable[i].t, aa, {ances, next_tip} };
            nodes[next_tip] = { death, next_node, {-1,-1}, label };
            nodes[ances].ances = next_node;
            const auto desc = nodes[aa].desc_idx(ances);
            nodes[aa].desc[desc] = next_node++;
            node_map[i] = next_tip++;
          }
        }
      }
      return tree_t{ ltable.age(), tips, true, std::move(nodes) };
    }


    tree_t to_tree(const ltable_t& ltable, std::false_type /*drop_extinct*/) {
      auto nodes = std::vector<node_t>(2 * ltable.size() - 1);
      nodes[ltable.size()] = { ltable[0].t, -1, {0, 1}, 0 };    // a.k.a. root
      nodes[0] = { ltable[0].death, ltable.size(), {-1,-1}, ltable[0].label };
      nodes[1] = { ltable[1].death, ltable.size(), {-1,-1}, ltable[1].label };
      int next_node = ltable.size() + 1;
      bool ultrametric = (ltable[0].death == 0.0) && (ltable[1].death == 0.0);
      for (int i = 2; i < ltable.size(); ++i) {
        const auto& e = ltable[i];
        const auto ances = e.ancestor;
        const auto aa = nodes[ances].ances;
        nodes[next_node] = { e.t, aa, { ances, i } };
        nodes[i] = { e.death, next_node, {-1, -1}, e.label };
        ultrametric &= (e.death == 0.0);
        nodes[ances].ances = next_node;
        const auto desc = nodes[aa].desc_idx(ances);
        nodes[aa].desc[desc] = next_node++;
      }
      return { ltable.age(), ltable.size(), ultrametric, std::move(nodes) };
    }

  }

  
  tree_t tree_t::from(const ltable_t& ltable, bool drop_extinct) {    
    return drop_extinct
      ? detail::to_tree(ltable, std::true_type{})
      : detail::to_tree(ltable, std::false_type{});
  }


  tree_t tree_t::from(const phylo_t& phylo) {
    auto e0 = phylo.e0();
    auto e1 = phylo.e1();
    auto len = phylo.len();
    auto tip_label = phylo.tip_label();
    auto nodes = std::vector<node_t>(e0.size() + 1);  // including root
    const auto tips = static_cast<int>(tip_label.size());
    double age = 0.0;
    nodes[tips].desc = {-1,-1};   // a.k.a. root
    for (int i = 0; i < static_cast<int>(e0.size()); ++i) {
      // cladewise traversal
      auto nn = e1[i] - 1;
      auto ances = e0[i] - 1;
      nodes[nn] = { nodes[ances].t + len[i], ances, {-1, -1} };
      auto desc = static_cast<int>(nodes[ances].desc[0] != -1);   // choose first non-assigned;
      nodes[ances].desc[desc] = nn;
      age = std::max(age, nodes[nn].t);
    }
    // some rounding errors might have creped in: rectify extremal times
    bool ultrametric = true;
    for (int i = 0; i < tips; ++i) {
      auto& node = nodes[i];
      if (std::abs(node.t - age) < present_tol) {
        node.t = age;
      }
      node.t = age - node.t;
      node.label = tip_label[i];
      ultrametric &= node.is_leaf() && (node.t == age);
    }
    for (int i = tips; i < static_cast<int>(nodes.size()); ++i) {
      nodes[i].t = age - nodes[i].t;
    }
    return tree_t{ age, static_cast<int>(tip_label.size()), ultrametric, std::move(nodes) };
  }

}
