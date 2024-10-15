/**
 * @file su2.hh
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief class representing an SU(2) matrix in the fundamental representation
 * @version 0.1
 * @date 2023-02-24
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <complex>
// #include<cmath>
#include <iostream>


#include "parameters.hh"
#include "accum_type.hh"
#include "dagger.hh"
#include "traceless_antiherm.hh"

using Complex = std::complex<double>;

/**
 * @brief representation of an SU(2) matrix in the fundamental representation with the
 * convention of eq. (4.23) of Gattringer&Lang
 * https://link.springer.com/book/10.1007/978-3-642-01850-3
 *
 */
class _su2 {
public:
  const size_t N_c = 2;
  
  
  explicit _su2() : a(0), b(0) {}
  explicit _su2(Complex a, Complex b) : a(a), b(b) {}
  _su2(const _su2 &U) : a(U.a), b(U.b) {}

  friend inline _su2 operator+(const _su2 &U1, const _su2 &U2);
  friend inline _su2 operator-(const _su2 &U1, const _su2 &U2);
  friend inline _su2 operator*(const _su2 &U1, const _su2 &U2);
  friend inline _su2 operator*(const Complex &U1, const _su2 &U2);
  friend inline _su2 operator*(const _su2 &U1, const Complex &U2);
  _su2 &operator*=(const _su2 &U1) {
    Complex a = this->a;
    this->a = a * U1.a - this->b * std::conj(U1.b);
    this->b = a * U1.b + this->b * std::conj(U1.a);
    return *this;
  }
  _su2 round(size_t n) const {
    double dn = n;
    return _su2(
      Complex(std::round(std::real(a) * dn), std::round(std::imag(a) * dn)) / dn,
      Complex(std::round(std::real(b) * dn), std::round(std::imag(b) * dn)) / dn);
  }

  inline Complex geta() const { return (a); }
  inline Complex getb() const { return (b); }
  inline void operator=(const _su2 &U) {
    a = U.geta();
    b = U.getb();
  }
  inline void operator+=(const _su2 &U) {
    a += U.a;
    b += U.b;
  }
  void set(const Complex _a, const Complex _b) {
    a = _a;
    b = _b;
  }
  void set_to_identity() {
    a = 1.0;
    b = 0.0;
  }
  inline _su2 dagger() const { return (_su2(std::conj(a), -b)); }
  inline double retrace() { return (2. * std::real(a)); }
  Complex det() { return (a * std::conj(a) + b * std::conj(b)); }
  void restoreSU() {
    double r = sqrt(std::abs(a) * std::abs(a) + std::abs(b) * std::abs(b));
    a /= r;
    b /= r;
  }

  void print() {
    std::cout << "--------------------\n";
    std::cout << a << " " << b << "\n";
    std::cout << -std::conj(b) << " " << std::conj(a) << "\n";
    std::cout << "--------------------\n";
  }

private:
  Complex a, b;
};

inline double retrace(_su2 const &U) {
  double a = std::real(U.geta());
  return (2 * a);
}

inline Complex trace(_su2 const &U) {
  double a = std::real(U.geta());
  return (Complex(2 * a, 0.));
}

template <> inline _su2 dagger(const _su2 &u) {
  const Complex a = u.geta();
  const Complex b = u.getb();
  _su2 udag(std::conj(a), -b);
  return udag;
}

template <> inline _su2 traceless_antiherm(const _su2 &x) {
  const std::complex<double> a1 = 0.5 * (x.geta() - std::conj(x.geta()));
  const std::complex<double> b1 = 0.5 * (x.getb() - std::conj(x.getb()));
  return (_su2(a1, b1));
}

inline _su2 operator*(const _su2 &U1, const _su2 &U2) {
  _su2 res;
  res.a = U1.a * U2.a - U1.b * std::conj(U2.b);
  res.b = U1.a * U2.b + U1.b * std::conj(U2.a);
  return (res);
}

inline _su2 operator+(const _su2 &U1, const _su2 &U2) {
  _su2 res;
  res.a = U1.a + U2.a;
  res.b = U1.b + U2.b;
  return (res);
}

inline _su2 operator-(const _su2 &U1, const _su2 &U2) {
  _su2 res;
  res.a = U1.a - U2.a;
  res.b = U1.b - U2.b;
  return (res);
}

inline _su2 operator*(const Complex &U1, const _su2 &U2) {
  _su2 res;
  res.a = U2.a * U1;
  res.b = U2.b * U1;
  return (res);
}

inline _su2 operator*(const _su2 &U1, const Complex &U2) {
  _su2 res;
  res.a = U1.a * U2;
  res.b = U1.b * U2;
  return (res);
}

using su2_accum = _su2;

inline _su2 accum_to_Group(const su2_accum &x) {
  _su2 U(x.geta(), x.getb());
  U.restoreSU();
  return U;
}

/**
 * @brief accumulation type for SU(2) matrices
 *
 * Incidentally, for SU(2) linear combinations of products its matrices can be still
 * parametrized by 2 numbers only. Therefore we use the same class for SU(2) elements and
 * their accumulations.
 *
 * @tparam
 */

template <> struct accum_type<_su2> {
  typedef _su2 type;
};

using su2 = _su2;
