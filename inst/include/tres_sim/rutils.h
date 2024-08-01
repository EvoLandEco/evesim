//  SPDX-FileCopyrightText: 2024 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
//  SPDX-License-Identifier: MIT

#ifndef RUTILS_H_INCLUDED
#define RUTILS_H_INCLUDED
#pragma once


// define to '0' to disable XPtr support
#ifndef RUTILS_TAGGED_XPTR
# define RUTILS_TAGGED_XPTR 1
#endif


#include "Rcpp.h"
#include "RcppParallel.h"
#include <cstdlib>
#include <random>
#include <string>
#include <memory>
#include <thread>


namespace rutils {


  // default random number generator
  using reng_t = std::mt19937_64;


  // creates an random number generator of type RENG.
  // more or less properly seeded by R's random number generator.
  // not ideal but ok for our purpose.
  template <typename RENG = reng_t>
  inline auto make_reng() {
    auto rseq = Rcpp::runif(5, 0, std::numeric_limits<uint32_t>::max());
    auto seq = std::seed_seq{
      static_cast<uint_least32_t>(rseq[0]),
      static_cast<uint_least32_t>(rseq[1]),
      static_cast<uint_least32_t>(rseq[2]),
      static_cast<uint_least32_t>(rseq[3]),
      static_cast<uint_least32_t>(rseq[4])
    };
    return RENG{seq};
  }


  // clean way to retrieve RcppParallel's concurrency setting
  // from RcppParallel::setThreadOptions(numThreads)
  inline int get_rcpp_num_threads() {
    auto* nt_env = std::getenv("RCPP_PARALLEL_NUM_THREADS");
    auto hc = static_cast<int>(std::thread::hardware_concurrency());
    int num_threads = (nullptr == nt_env) ? hc : std::atoi(nt_env);
    return std::clamp(num_threads, 1, hc);
  }


  // RAII scoped guard for tbb::global_control
  class tbb_global_control_guard {
  public:
    // creates tbb::global_control iff 
    //  1 < get_rcpp_num_threads()
    //  1 == get_rcpp_num_threads() && true == single_threaded_control_guard
    explicit tbb_global_control_guard(bool single_threaded_control_guard = true)
      : num_threads_ (get_rcpp_num_threads()) {
      if ((num_threads_ != 1) || single_threaded_control_guard) {
        gc_.reset(new tbb::global_control(tbb::global_control::max_allowed_parallelism, num_threads_));
      }
    }

    int num_threads() const noexcept { return num_threads_; }

  private:
    int num_threads_;
    std::unique_ptr<tbb::global_control> gc_;
  };


#if RUTILS_TAGGED_XPTR == 1


  template <typename T>
  inline bool is_tagged_xptr(SEXP x, const char* Tag) {
    if (TYPEOF(x) == EXTPTRSXP) {
      if (auto tag = R_ExternalPtrTag(x); tag != nullptr) {
        if (Rcpp::is<Rcpp::String>(tag)) {
          auto rtag = Rcpp::as<Rcpp::String>(tag);
          return (rtag == Tag);
        }
      }
    }
    return false;
  }


  // returns a copy of the tagged XPtr<T>
  template <typename T>
  inline Rcpp::XPtr<T> tagged_xptr(SEXP x) {
    if (is_tagged_xptr<T>(x)) {
      return Rcpp::XPtr<T>(x, R_ExternalPtrTag(x), R_NilValue);
    }
    throw std::runtime_error("external pointer tag mismatch");
    return Rcpp::XPtr<T>(R_NilValue);
  }


  // returns a copy of the tagged XPtr<T>
  template <typename T>
  inline Rcpp::XPtr<T> tagged_xptr(SEXP x, const char* Tag) {
    if (is_tagged_xptr<T>(x, Tag)) {
      return Rcpp::XPtr<T>(x, R_ExternalPtrTag(x), R_NilValue);
    }
    throw std::runtime_error(std::string("external pointer mismatch: expected ") + Tag);
    return Rcpp::XPtr<T>(R_NilValue);
  }


  // creates a tagged XPtr<T>
  template <typename T>
  inline Rcpp::XPtr<T> tagged_xptr(T* ptr, const char* Tag) {
    return Rcpp::XPtr<T>(ptr, true, Rcpp::wrap(Tag));
  }


#endif  // RUTILS_TAGGED_XPTR == 1


  template <typename T>
  class RView {
  public:
    RView(T* first, size_t n, size_t m = 1) : begin_(first), n_(n), m_(m) {}
    size_t size() const noexcept { return n_ * m_; }
    size_t length() const noexcept { return n_ * m_; }
    size_t nrow() const noexcept { return n_; }
    size_t ncol() const noexcept { return m_; }

    T* begin() noexcept { return begin_; }
    T* end() noexcept { return begin_ + size(); }

    // implicit conversions
    operator RcppParallel::RVector<T>() const { return RcppParallel::RVector<T>(*this); }
    operator RcppParallel::RMatrix<T>() const { return RcppParallel::RMatrix<T>(*this); }
    
  private:
    T* begin_;
    size_t n_, m_;
  };


} // namespace rutils

#endif
