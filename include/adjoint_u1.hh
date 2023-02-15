/**
 * @file adjoint_u1.hh
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
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

template <typename Float = double> inline adjointu1<Float> get_deriv(Complex &A) {
  return adjointu1<Float>(std::imag(A));
}
