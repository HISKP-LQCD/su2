/**
 * @file adjoint_su3.hh
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief class and routines for an su(3) matrix in the adjoint representation
 * @version 0.1
 * @date 2023-02-15
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include "geometry.hh"
#include "su3.hh"

#include <array>
#include <cassert>
#include <cmath>
#include <random>
#include <vector>

/**
 * @brief element of the su(3) algebra
 *
 * the element is parametrized by 8 real numbers (number of generators for SU(3))
 *
 * @tparam Float
 */
template <typename Float> class adjointsu3 {
public:
  adjointsu3(const std::array<Float, 8> &_arr) : arr(_arr) {}
  adjointsu3() { this->setzero(); }
  void flipsign() {
    for (size_t i = 0; i < 8; i++) {
      arr[i] = -arr[i];
    }
  }
  std::array<Float, 8> get_arr() const { return arr; }
  void set_arr(const std::array<Float, 8> &_arr) { arr = _arr; }
  void setzero() {
    const std::array<Float, 8> _arr{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    arr = _arr;
  }
  adjointsu3<Float> round(const size_t &n) const {
    Float dn = n;
    std::array<Float, 8> arr2;
    for (size_t i = 0; i < 8; i++) {
      arr2[i] = std::round(arr[i] * dn) / dn;
    }
    return adjointsu3(arr2);
  }

  void operator=(const adjointsu3 &A) { arr = A.get_arr(); }
  void operator+=(const adjointsu3 &A) {
    const std::array<Float, 8> arr_A = A.get_arr();
    for (size_t i = 0; i < 8; i++) {
      arr[i] += arr_A[i];
    }
  }
  void operator-=(const adjointsu3 &A) {
    const std::array<Float, 8> arr_A = A.get_arr();
    for (size_t i = 0; i < 8; i++) {
      arr[i] -= arr_A[i];
    }
  }

private:
  std::array<Float, 8> arr;
};

template <typename Float = double> inline adjointsu3<Float> get_deriv(const su3 &A) {
  const std::array<std::complex<Float>, 9> arr_A = A.get_arr();
  const std::array<Float, 8> arr;
  for (size_t i = 0; i < 8; i++) {
    arr[i] = 2.0 * std::imag(arr_A[i+1]);
  }
  return adjointsu3<Float>(arr);
}
