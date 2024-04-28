/**
 * @file adjoint_u1.hh
 * @author Simone Romiti (simone.romiti.1994@gmail.com)
 * @brief class and routines for an u(1) matrix in the adjoint representation
 * (constant*Identity)
 * @version 0.1
 * @date 2023-02-15
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include "geometry.hh"
#include "u1.hh"

#include <array>
#include <cassert>
#include <cmath>
#include <random>
#include <vector>

/**
 * @brief element of the adjount representation of U(1), i.e. just a number proportional
 * to the Identity.
 *
 * Note: We use the convention such that tr(t_i*t_j) = delta_{ij} (no factor 1/2), so the
 * generator of the algebra is the Identity, not divided by sqrt(2)
 *
 * @tparam Float
 */
template <typename Float> class adjointu1 {
public:
  adjointu1(Float _a) : a(_a) {}
  adjointu1() : a(0.) {}
  void flipsign() { a = -a; }
  Float geta() const { return a; }
  void seta(Float _a) { a = _a; }
  void setzero() { a = 0.; }

  adjointu1<Float> round(size_t n) const {
    Float dn = n;
    return adjointu1(std::round(a * dn) / dn);
  }
  void operator=(const adjointu1 &A) { a = A.geta(); }
  void operator+=(const adjointu1 &A) { a += A.geta(); }
  void operator-=(const adjointu1 &A) { a -= A.geta(); }

private:
  Float a;
};

template <typename Float = double> inline adjointu1<Float> get_deriv(const u1_accum &A) {
  return adjointu1<Float>(std::imag(A));
}

template <typename Float = double> inline Float get_abs(const u1_accum &A) {
  return Float(std::abs(A));
}

template <typename Float = double> inline Float get_phase(const u1_accum &A) {
  return Float(std::arg(A));
}

template <typename Float>
inline adjointu1<Float> operator*(const Float &x, const adjointu1<Float> &A) {
  return adjointu1<Float>(x * A.geta());
}
