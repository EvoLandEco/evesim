//  SPDX-FileCopyrightText: 2024 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
//  SPDX-License-Identifier: MIT

#ifndef TREE_VISIT_H_INCLUDED
#define TREE_VISIT_H_INCLUDED
#pragma once

#include <stack>
#include <queue>
#include "tres_sim/tree.h"


namespace tres_sim {


   // traversal types for 'visit'
  namespace visit_order {
    struct depth_first {};    // tag depth first, pre-order traversal
    struct breadth_first {};  // tag breadth-first traversal
    struct clade_wise {};     // tag clade-wise (depth-first per crown clade)
    struct linear {};         // tag linear traversal (iterating over valid nodes)
  }


  // calls 'fun(idx)' for each node visited while traversing
  // the tree in depth-first order, starting from 'start_idx' 
  // until 'fun' returns 'true'.
  // Note that this function accepts Fun->void as well as Fun->bool
  template <typename Fun>
  inline bool visit(visit_order::depth_first, const std::vector<node_t>& nodes, int start_index, Fun&& fun) {
    if (nodes.empty()) return false;
    auto stack = std::vector<int>{};
    stack.push_back(start_index);
    while (!stack.empty()) {
      auto idx = stack.back();
      stack.pop_back();
      if constexpr (std::is_same_v<bool , std::invoke_result_t<Fun, int>>) {
        if (fun(idx)) return true;
      }
      else {
        fun(idx);
      }
      if (nodes[idx].is_inode()) {
        stack.push_back(nodes[idx].desc[1]);
        stack.push_back(nodes[idx].desc[0]);
      }
    }
    return false;
  }


  // calls 'fun(idx)' for each node visited while traversing
  // the tree in depth-first order, starting from 'root' 
  // until 'fun' returns 'true'.
  // Note that this function accepts Fun->void as well as Fun->bool
  template <typename Fun>
  inline bool visit(visit_order::depth_first, const tree_t& tree, Fun&& fun) {
    return visit(visit_order::depth_first{}, tree.nodes, tree.root, std::forward<Fun>(fun));
  }


  // calls 'fun(idx)' for each node visited while traversing
  // the tree in breadth-first order, starting from 'start_idx' 
  // until 'fun' returns 'true'.
  // Note that this function accepts Fun->void as well as Fun->bool
  template <typename Fun>
  inline bool visit(visit_order::breadth_first, const std::vector<node_t>& nodes, int start_index, Fun&& fun) {
    if (nodes.empty()) return false;
    auto queue = std::queue<int>{};
    queue.push(start_index);
    while (!queue.empty()) {
      auto idx = queue.front();
      queue.pop();
      if constexpr (std::is_same_v<bool , std::invoke_result_t<Fun, int>>) {
        if (fun(idx)) return true;
      }
      else {
        fun(idx);
      }
      if (nodes[idx].is_inode()) {
        queue.push(nodes[idx].desc[1]);
        queue.push(nodes[idx].desc[0]);
      }
    }
    return false;
  }


  // calls 'fun(idx)' for each node visited while traversing
  // the tree in breadth-first order, starting from 'root' 
  // until 'fun' returns 'true'.
  // Note that this function accepts Fun->void as well as Fun->bool
  template <typename Fun>
  inline bool visit(visit_order::breadth_first, const tree_t& tree, Fun&& fun) {
    return visit(visit_order::breadth_first{}, tree.nodes, tree.root, std::forward<Fun>(fun));
  }


  // calls 'fun(idx)' for each node visited while traversing
  // the tree in clade-first order, until 'fun' returns 'true'.
  // Note that this function accepts Fun->void as well as Fun->bool
  template <typename Fun>
  inline bool visit(visit_order::clade_wise, const tree_t& tree, Fun fun) {
    auto& nroot = tree.nodes[tree.root];
    return visit(visit_order::depth_first{}, tree.nodes, nroot.desc[0], fun) ||
           visit(visit_order::depth_first{}, tree.nodes, nroot.desc[1], fun);
  }


  // calls 'fun(idx)' for each node visited while traversing
  // the tree in linear oder until 'fun' returns 'true'.
  // Note that this function accepts Fun->void as well as Fun->bool
  template <typename Fun>
  inline bool visit(visit_order::linear, const std::vector<node_t>& nodes, int start_index, Fun&& fun) {
    for (int idx = start_index; idx < static_cast<int>(nodes.size()); ++idx) {
      if constexpr (std::is_same_v<bool , std::invoke_result_t<Fun, int>>) {
        if (fun(idx)) return true;
      }
      else {
        fun(idx);
      }
    }
    return false;
  }


  // calls 'fun(idx)' for each node visited while traversing
  // the tree in depth-first order, starting from 'root' 
  // until 'fun' returns 'true'.
  // Note that this function accepts Fun->void as well as Fun->bool
  template <typename Fun>
  inline bool visit(visit_order::linear, const tree_t& tree, Fun&& fun) {
    return visit(visit_order::linear{}, tree.nodes, tree.root, std::forward<Fun>(fun));
  }


}
#endif
