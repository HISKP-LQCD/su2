/**
 * @file su3.hh
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief class and routines for an SU(3) matrix in the fundamental representation
 * @version 0.1
 * @date 2023-02-15
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <complex>
// #include<cmath>
#include <iostream>

#include "accum_type.hh"
#include "dagger.hh"
#include "traceless_antiherm.hh"

using Complex = std::complex<double>;

/**
 * @brief representation of an SU(3) matrix as U = a_0*1 + \sum_{i=1}^{8} a_i *
 \lambda_i, where lambda_i are the Gell-Mann matrices normalized such that
 tr(lambda_i*lambda_j)=2*delta_ij
 *
 */
class _su3 {
public:
  const size_t N_c = 3;
  explicit _su3() {
    const std::array<Complex, 9> _arr{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    arr = _arr;
  }
  explicit _su3(const std::array<Complex, 9> &_arr) : arr(_arr) {}
  _su3(const _su3 &U) : arr(U.get_arr()) {}

  // friend inline _su3 operator+(const _su3 &U1, const _su3 &U2);
  // friend inline _su3 operator-(const _su3 &U1, const _su3 &U2);
  friend inline _su3 operator*(const _su3 &U1, const _su3 &U2);
  friend inline _su3 operator*(const Complex &U1, const _su3 &U2);
  friend inline _su3 operator*(const _su3 &U1, const Complex &U2);
  // _su3 &operator*=(const _su3 &U1) {
  //   // CHANGE
  //   // Complex a = this->a;
  //   // this->a = a * U1.a - this->b * std::conj(U1.b);
  //   // this->b = a * U1.b + this->b * std::conj(U1.a);
  //   // return *this;
  // }
  _su3 round(size_t n) const {
    double dn = n;
    std::array<Complex, 9> arr2;
    for (size_t i = 0; i < 9; i++) {
      arr2[i] = Complex(std::round(std::real(arr[i]) * dn) / dn,
                        std::round(std::imag(arr[i]) * dn) / dn);
    }
    return _su3(arr2);
  }

  inline std::array<Complex, 9> get_arr() const { return arr; }
  inline void operator=(const _su3 &U) { arr = U.get_arr(); }
  // inline void operator+=(const _su3 &U) {
  //   a += U.a;
  //   b += U.b;
  // }
  void set(const std::array<Complex, 9> &_arr) { arr = _arr; }
  inline _su3 dagger() const {
    std::array<Complex, 9> arr2;
    for (size_t i = 0; i < 9; i++) {
      arr2[i] = std::conj(arr[i]);
    }
    return _su3(arr2);
  }
  inline Complex trace() { return double(N_c) * arr[0]; }
  inline double retrace() { return std::real(this->trace()); }
  // Complex det() { return (a * std::conj(a) + b * std::conj(b)); }
  // void restoreSU() {
  //   double r = sqrt(std::abs(a) * std::abs(a) + std::abs(b) * std::abs(b));
  //   a /= r;
  //   b /= r;
  // }

private:
  std::array<Complex, 9> arr; // {a_0, a_1, ..., a_8}
};

inline Complex trace(_su3 const &U) {
  return double(U.N_c) * U.get_arr()[0];
}

inline double retrace(_su3 const &U) {
  return std::real(trace(U));
}

template <> inline _su3 dagger(const _su3 &u) {
  return u.dagger();
}

template <> inline _su3 traceless_antiherm(const _su3 &x) {
  std::array<Complex, 9> arr = x.get_arr();
  // make it anti-hermitian (loop can start from 1)
  for (size_t i = 1; i < 9; i++) {
    arr[i] = (arr[i] - std::conj(arr[i]));
  }
  arr[0] = 0.0; // make it traceless (trace comes from the 0-th element only)
  return _su3(arr);
}

// inline _su3 operator*(const _su3 &U1, const _su3 &U2) {
//   _su3 res;
//   res.a = U1.a * U2.a - U1.b * std::conj(U2.b);
//   res.b = U1.a * U2.b + U1.b * std::conj(U2.a);
//   return (res);
// }

// inline _su3 operator+(const _su3 &U1, const _su3 &U2) {
//   _su3 res;
//   res.a = U1.a + U2.a;
//   res.b = U1.b + U2.b;
//   return (res);
// }

// inline _su3 operator-(const _su3 &U1, const _su3 &U2) {
//   _su3 res;
//   res.a = U1.a - U2.a;
//   res.b = U1.b - U2.b;
//   return (res);
// }

// inline _su3 operator*(const Complex &U1, const _su3 &U2) {
//   _su3 res;
//   res.a = U2.a * U1;
//   res.b = U2.b * U1;
//   return (res);
// }

// inline _su3 operator*(const _su3 &U1, const Complex &U2) {
//   _su3 res;
//   res.a = U1.a * U2;
//   res.b = U1.b * U2;
//   return (res);
// }

using su3 = _su3;
