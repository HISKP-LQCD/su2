/**
 * @file adjoint_su2.hh
 * @author Simone Romiti (simone.romiti.1994@gmail.com)
 * @brief class and routines for an su(2) matrix in the adjoint representation
 * @version 0.1
 * @date 2023-02-15
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include "geometry.hh"
#include "su2.hh"

#include <array>
#include <cassert>
#include <cmath>
#include <random>
#include <vector>

template <typename Float> class adjointsu2 {
public:
  adjointsu2(Float _a, Float _b, Float _c) : a(_a), b(_b), c(_c) {}
  adjointsu2() : a(0.), b(0.), c(0.) {}
  void flipsign() {
    a = -a;
    b = -b;
    c = -c;
  }
  Float geta() const { return a; }
  Float getb() const { return b; }
  Float getc() const { return c; }
  void seta(Float _a) { a = _a; }
  void setb(Float _a) { b = _a; }
  void setc(Float _a) { c = _a; }
  void setzero() { a = b = c = 0.; }
  adjointsu2<Float> round(size_t n) const {
    Float dn = n;
    return adjointsu2(std::round(a * dn) / dn, std::round(b * dn) / dn,
                      std::round(c * dn) / dn);
  }
  void operator=(const adjointsu2 &A) {
    a = A.geta();
    b = A.getb();
    c = A.getc();
  }
  void operator+=(const adjointsu2 &A) {
    a += A.geta();
    b += A.getb();
    c += A.getc();
  }
  void operator-=(const adjointsu2 &A) {
    a -= A.geta();
    b -= A.getb();
    c -= A.getc();
  }

private:
  Float a, b, c;
};

template <typename Float = double>
inline adjointsu2<Float> get_deriv(const su2_accum &A) {
  const Complex a = A.geta(), b = A.getb();
  return adjointsu2<Float>(2. * std::imag(b), 2. * std::real(b), 2. * std::imag(a));
}

template <typename Float>
inline adjointsu2<Float> operator*(const Float &x, const adjointsu2<Float> &A) {
  return adjointsu2<Float>(x * A.geta(), x * A.getb(), x * A.getc());
}
