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
 * any su(3) matrix can be expressed as the anti-hermitian traceless part of an SU(3)
 * matrix therefore, the object is parametrized by 8 real numbers (number of generators
 * for SU(3))
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

/**
 * @brief returns (-i)*(B - B.dagger()), where B = A - (tr(a)/N_c) * 1, parametrized
 * as an element of su(3)
 *
 * @tparam Float
 * @param A
 * @return adjointsu3<Float>
 */
template <typename Float = double> inline adjointsu3<Float> get_deriv(const su3 &A) {
//  std::cout << "ERROR: please implement " << __func__ << " for SU(3) \n";
//  std::abort();
  const su3 A_thh = traceless_antiherm(A);
  const std::array<Complex, 3> u = A_thh.get_u();
  const std::array<Complex, 3> v = A_thh.get_v();
  const std::array<Float, 8> arr = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // = ???
  return adjointsu3<Float>(arr);
}

template <typename Float>
inline adjointsu3<Float> operator*(const Float &x, const adjointsu3<Float> &A) {
  std::array<Float, 8> arr = A.get_arr();
  for (size_t i = 0; i < 8; i++) {
    arr[i] *= x;
  }
  return adjointsu3(arr);
}
