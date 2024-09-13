//  SPDX-FileCopyrightText: 2024 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
//  SPDX-License-Identifier: MIT

#ifndef EVE_SIM_H_INCLUDED
#define EVE_SIM_H_INCLUDED
#pragma once

#include <algorithm>
#include <numeric>
#include <functional>
#include <tuple>
#include "evesim/phylo.h"


namespace eve { 

/*   template <typename EVENTGEN>
  inline Ltable dd_sim(EVENTGEN& event_gen, double age, size_t retry);


  struct Param {
    explicit Param(const Rcpp::NumericVector& rparam) noexcept
    : la(rparam[0]), mu(rparam[1]), K(rparam[2]), r((rparam.size() < 4) ? 0.0 : rparam[3]) 
    {}

    double la;
    double mu;
    double K;
    double r;
  };

  struct lamu { 
    double laN;
    double muN;
  };

  using lamuN_fun = std::function<lamu(double)>;

  inline lamuN_fun lamuN(const Param& p, double ddmodel) 
  {
    const auto imodel = static_cast<int>(100 * ddmodel);  // leave some room for future sub-models.
    switch (imodel) {
      case 100: return [=](double N) { return lamu{ std::max(0.0, p.la - (p.la - p.mu) * N / p.K), p.mu }; }; break;
      case 130: return [=](double N) { return lamu{ std::max(0.0, p.la * (1.0 - N / p.K)), p.mu }; }; break;
      case 200:
      case 210:
      case 220: {
        double al=(std::log(p.la/p.mu)/std::log(p.K + 1.0));
        return [=](double N) { return lamu{ std::pow(p.la * (N + 1.0), -al), p.mu }; }; 
        break;
      }
      case 230: return [=](double N) { return lamu{ std::pow(p.la * (N + 1.0), -p.K), p.mu }; }; break;
      case 300: return [=](double N) { return lamu{ p.la, p.mu + (p.la + p.mu) * N / p.K }; }; break;
      case 400:
      case 410: {
        double al=(std::log(p.la / p.mu) / std::log(p.K + 1.0));
        return [=](double N) { return lamu{ p.la, std::pow(p.mu * (N + 1.0), al) }; };
        break; 
      }
      case 420: return [=](double N) { return lamu{ p.la, p.mu }; }; break;
      case 500: {
        double t0 = p.la - 1.0 / (p.r + 1.0) * (p.la - p.mu) / p.K;
        double t1 = p.r / (p.r + 1.0) * (p.la - p.mu) / p.K;
        return [=](double N) { return lamu{ std::max(0.0, t0 * N), p.mu + t1 * N }; }; break;
      }
      default: throw std::runtime_error("unknown ddmodel");
    }
  }


  // event generator for standard dd_sim
  class lamuN_event_gen {
  public:
    lamuN_event_gen(lamuN_fun&& lamuN) : 
      lamuN_(std::move(lamuN)),
      reng_(make_reng<reng_t>())
    {}

    struct event_t {
      bool is_speciation;
      double t;             // time
      int ranI;             // the affected species
    };

    event_t operator()(double t, size_t tips) {
      const auto [laN, muN] = lamuN_(tips);
      return { 
        std::bernoulli_distribution{laN / (laN + muN)}(reng_),                               // true iff speciation
        t - std::exponential_distribution{(laN + muN) * static_cast<double>(tips)}(reng_),   // time
        std::uniform_int_distribution<>(0, static_cast<int>(tips - 1))(reng_)                // ranI, the affected species
      }; 
    }
  
    void reset () const noexcept {}

  private:
    lamuN_fun lamuN_;
    reng_t reng_;
  };


  template <typename EVENTGEN>
  inline Ltable dd_sim(EVENTGEN& event_gen, double age, size_t retry) 
  {
    auto LL = Ltable{};
    auto tips = std::vector<int>{};   // indices of extant species
    for (; retry; --retry) {
      LL.assign({
        Ltable_entry{age, -1, -1.0, true}, 
        Ltable_entry{age, -1, -1.0, false}
      });
      tips.assign({0, 1});   // list of extant species
      event_gen.reset();
      auto event = event_gen(age, tips.size());
      std::array<int, 2> ncl{1, 1};    // # species per crown lineage
      auto cl_present = [&]() { return 0 != (ncl[0] * ncl[1]); };  // descendants from both crown lineages present?
      while (cl_present() && (event.t >= 0.0)) {
        const auto ranI = event.ranI;
        const auto ranL = tips[ranI];
        const auto sranL = LL[ranL].flag;    // inherited crown lineage flag of ancestor
        if (event.is_speciation) {
          // speciation
          ++ncl[sranL];
          tips.push_back(static_cast<int>(LL.size()));
          LL.push_back({ event.t, ranL, -1, sranL });
        }
        else {
          // extinction
          LL[ranL].death = event.t;
          --ncl[sranL];
          tips[ranI] = tips.back();
          tips.pop_back();
        }
        event = event_gen(event.t, tips.size());
      }
      if (cl_present()) {
        return LL;    // success.
      }
    }
    return Ltable{};  // failure: retry-overrun.
  }


  // Example of how to use a tree instead of a Ltable internally in simulations.
  // Callees expect Ltables, thus this is quite stupid for now.
  template <typename EVENTGEN>
  inline Ltable dd_sim_tree(EVENTGEN& event_gen, double age, size_t retry) 
  {
    struct fnode_t : node_t {
      int flag = 0;             // crown lineage flag
    };
    auto nodes = std::vector<fnode_t>{};
    auto flags = std::vector<int>{};
    auto tips = std::vector<int>{};    // iff t == -1.0;
    for (; retry; --retry) {
      nodes.assign({
        {age, -1, {1, 2}},            // root
        {-1.0, 0, {-1, -1}, 1},       // clade 0 (see README)
        {-1.0, 0, {-1, -1}, 0}        // clade 1
      });
      tips.assign({1,2});
      event_gen.reset();
      auto event = event_gen(age, tips.size());
      std::array<int, 2> ncl{1, 1};    // # species per crown lineage
      auto cl_present = [&]() { return 0 != (ncl[0] * ncl[1]); };  // descendants from both crown lineages present?
      while (cl_present() && (event.t >= 0.0)) {
        const auto leaf = tips[event.ranI];
        auto& node = nodes[leaf];
        node.t = event.t;
        const auto flag = node.flag;   // inherited
        if (event.is_speciation) {
          // speciation
          const auto nn = static_cast<int>(nodes.size());
          nodes.push_back({ -1.0, leaf, {-1, -1}, flag }); 
          nodes.push_back({ -1.0, leaf, {-1, -1}, flag }); 
          node.desc = {nn, nn + 1};
          tips[event.ranI] = nn;
          tips.push_back(nn + 1);
          ++ncl[flag];
        }
        else {
          // extinction
          --ncl[node.flag];
          tips[event.ranI] = tips.back();
          tips.pop_back();
        }
        event = event_gen(event.t, tips.size());
      }
      if (cl_present()) {
          auto tree = tree_t<fnode_t>{ .age = age, .root = 0, .nnode = static_cast<int>(nodes.size()) - (ncl[0] + ncl[1]), .nodes = std::move(nodes) };
          return to_Ltable(tree);
      }
    }
    return Ltable{};  // failure: retry-overrun.
  }
 */
} // namespce eve

#endif
