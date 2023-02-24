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
#include "su3_accum.hh"

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
 * @brief returns A = B - (tr(B)/N_c) * 1_{3x3}, where B = U - U^{\dagger}.
 * The matrix A is returned as an element of the algebra su(3), i.e. giving the
 * coefficients in front of the Gell-Mall matrices
 *
 * @tparam Float
 * @param A
 * @return adjointsu3<Float>
 */
template <class Float = double> inline adjointsu3<Float> get_deriv(const su3_accum &U) {
  su3_accum A = U - U.dagger();
  _su3 Id; // by default is the identity
  A = A - (trace(A) / double(U.N_c)) * Id;
  const std::array<Complex, 3> u = A.get_u();
  const std::array<Complex, 3> v = A.get_v();
  const std::array<Complex, 3> w = A.get_w();

  std::array<Float, 8> arr;

  arr[0] = std::imag(u[1]) + std::imag(v[0]);
  arr[1] = std::real(u[1]) - std::real(v[0]);
  arr[2] = std::imag(u[0]) - std::imag(v[1]);
  arr[3] = std::imag(u[2]) + std::imag(w[0]);
  arr[4] = std::real(u[2]) - std::real(w[0]);
  arr[5] = std::imag(v[2]) + std::imag(w[1]);
  arr[6] = std::real(v[2]) - std::real(w[1]);
  arr[7] = sqrt(3.0) * std::imag(u[0]) / 3.0 + sqrt(3.0) * std::imag(v[1]) / 3.0 -
           2.0 * sqrt(3.0) * std::imag(w[2]) / 3.0;

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
