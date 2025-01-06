/**
 * @file su3_accum.hh
 * @author Simone Romiti (simone.romiti.1994@gmail.com)
 * @brief accumulation type for SU(3) matrices
 * @version 0.1
 * @date 2023-02-24
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <array>
#include <complex>
#include <iostream>

#include <Eigen/Dense> // Include Eigen library

#include "accum_type.hh"
#include "dagger.hh"
#include "su3.hh"
#include "traceless_antiherm.hh"

/**
 * @brief class describing a 3x3 matrix with rows (u, v, w). Its used to accumulate sums
 * and differences of plaquettes and/or staples (which is not unitary anymore)
 *
 */
class _su3_accum {
public:
  const size_t N_c = 3;
  explicit _su3_accum() {
    this->set_to_zero(); // default: zero matrix
  }
  explicit _su3_accum(const std::array<Complex, 3> &_u,
                      const std::array<Complex, 3> &_v,
                      const std::array<Complex, 3> &_w)
    : u(_u), v(_v), w(_w) {}
  _su3_accum(const _su3_accum &U) : u(U.get_u()), v(U.get_v()), w(U.get_w()) {}
  _su3_accum(const _su3 &U) : u(U.get_u()), v(U.get_v()), w(U.get_w()) {}

  friend inline _su3_accum operator+(const _su3_accum &U1, const _su3_accum &U2);
  friend inline _su3_accum operator-(const _su3_accum &U1, const _su3_accum &U2);
  friend inline _su3 operator*(const _su3 &U1, const _su3 &U2);

  template <class matr_33> _su3_accum &operator*=(const matr_33 &U2) {
    const std::array<Complex, 3> u1 = (*this).u, v1 = (*this).v, w1 = (*this).w;
    const std::array<Complex, 3> u2 = U2.get_u(), v2 = U2.get_v(), w2 = U2.get_w();
    for (size_t i = 0; i < this->N_c; i++) {
      (*this).u[i] = u1[0] * u2[i] + u1[1] * v2[i] + u1[2] * w2[i];
      (*this).v[i] = v1[0] * u2[i] + v1[1] * v2[i] + v1[2] * w2[i];
      (*this).w[i] = w1[0] * u2[i] + w1[1] * v2[i] + w1[2] * w2[i];
    }
    return *this;
  }

  inline std::array<Complex, 3> get_u() const { return u; }
  inline std::array<Complex, 3> get_v() const { return v; }
  inline std::array<Complex, 3> get_w() const { return w; }

  inline void operator=(const _su3 &U) {
    u = U.get_u();
    v = U.get_v();
    w = U.get_w();
  }

  inline void operator=(const _su3_accum &U) {
    u = U.get_u();
    v = U.get_v();
    w = U.get_w();
  }

  inline void operator+=(const _su3_accum &U) {
    for (size_t i = 0; i < N_c; i++) {
      u[i] += U.u[i];
      v[i] += U.v[i];
      w[i] += U.w[i];
    }
  }
  inline void operator-=(const _su3_accum &U) {
    for (size_t i = 0; i < N_c; i++) {
      u[i] -= U.u[i];
      v[i] -= U.v[i];
      w[i] -= U.w[i];
    }
  }
  void set(const std::array<Complex, 3> &_u,
           const std::array<Complex, 3> &_v,
           const std::array<Complex, 3> &_w) {
    u = _u;
    v = _v;
    w = _w;
  }
  void set_to_zero() {
    u = {0.0, 0.0, 0.0};
    v = u;
    w = u;
  }

  inline _su3_accum dagger() const {
    const std::array<Complex, 3> u2 = {std::conj(u[0]), std::conj(v[0]), std::conj(w[0])},
                                 v2 = {std::conj(u[1]), std::conj(v[1]), std::conj(w[1])},
                                 w2 = {std::conj(u[2]), std::conj(v[2]), std::conj(w[2])};
    return _su3_accum(u2, v2, w2);
  }
  inline Complex trace() const { return (u[0] + v[1] + w[2]); }
  inline double retrace() const { return std::real(this->trace()); }
  inline Complex det() const {
    Complex res = 0.0;
    res += +u[0] * (v[1] * w[2] - v[2] * w[1]);
    res += -u[1] * (v[0] * w[2] - v[2] * w[0]);
    res += +u[2] * (v[0] * w[1] - v[1] * w[0]);
    return res;
  }

  // returns V and D, where M is the original matrix and M=V*D*V^dagger
  // NOTE: this method is well defined only if (*this) is a Hermitean matrix
  inline std::vector<_su3_accum> diagonalize() {
    // Convert the _su3_accum into an Eigen matrix
    Eigen::Matrix3cd matrix;
    matrix << u[0], u[1], u[2], v[0], v[1], v[2], w[0], w[1], w[2];

    // Perform Eigen decomposition for Hermitian matrix
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd> solver(matrix);

    // Extract eigenvalues (diagonal elements of D) and eigenvectors (columns of V)
    Eigen::Vector3d eigenvalues = solver.eigenvalues();
    Eigen::Matrix3cd eigenvectors = solver.eigenvectors();

    // Convert V and D back to _su3_accum
    const std::array<Complex, 3> V1 = {eigenvectors(0, 0), eigenvectors(1, 0),
                                       eigenvectors(2, 0)};
    const std::array<Complex, 3> V2 = {eigenvectors(0, 1), eigenvectors(1, 1),
                                       eigenvectors(2, 1)};
    const std::array<Complex, 3> V3 = {eigenvectors(0, 2), eigenvectors(1, 2),
                                       eigenvectors(2, 2)};
    _su3_accum V(V3, V2, V3);

    const std::array<Complex, 3> D1 = {eigenvalues(0), 0.0, 0.0};
    const std::array<Complex, 3> D2 = {0.0, eigenvalues(1), 0.0};
    const std::array<Complex, 3> D3 = {0.0, 0.0, eigenvalues(2)};

    _su3_accum D(D1, D2, D3);

    // Return a vector of _su3_accum: [V, D]
    return {V, D};
  }

  // return the matrix with the diagonal elements D_{ii} raised to the power alpha
  inline _su3_accum pow_diagonal(const double &alpha) const {
    _su3_accum D_alpha = (*this);
    D_alpha.set({std::pow(u[0], alpha), 0.0, 0.0}, {0.0, std::pow(v[1], alpha), 0.0},
                {0.0, 0.0, std::pow(w[2], alpha)});
    return D_alpha;
  }

  inline su3 to_SU3() const {
    const su3 U(u, v);
    return U;
  }

  void print() const {
    std::cout << "----------------------------------\n";
    std::cout << u[0] << " " << u[1] << " " << u[2] << "\n";
    std::cout << v[0] << " " << v[1] << " " << v[2] << "\n";
    std::cout << w[0] << " " << w[1] << " " << w[2] << "\n";
    std::cout << "----------------------------------\n";
  }

private:
  std::array<Complex, 3> u, v, w;
};

inline _su3_accum operator+(const _su3_accum &U1, const _su3_accum &U2) {
  _su3_accum U = U1;
  U += U2;
  return U;
}

inline _su3_accum operator-(const _su3_accum &U1, const _su3_accum &U2) {
  _su3_accum U = U1;
  U -= U2;
  return U;
}

inline Complex trace(_su3_accum const &U) {
  return U.trace();
}

inline double retrace(_su3_accum const &U) {
  return U.retrace();
}

template <> inline _su3_accum dagger(const _su3_accum &u) {
  return u.dagger();
}

inline _su3_accum operator*(const _su3_accum &U1, const _su3_accum &U2) {
  _su3_accum U = U1;
  U *= U2;
  return U;
}

inline _su3_accum operator*(const double &z, const _su3_accum &U) {
  std::array<Complex, 3> u = U.get_u(), v = U.get_v(), w = U.get_w();
  for (size_t i = 0; i < U.N_c; i++) {
    u[i] *= double(z);
    v[i] *= double(z);
    w[i] *= double(z);
  }
  return _su3_accum(u, v, w);
}

inline _su3_accum operator*(const _su3_accum &U, const double &z) {
  return z * U;
}

inline _su3_accum operator*(const Complex &z, const _su3_accum &U) {
  std::array<Complex, 3> u = U.get_u(), v = U.get_v(), w = U.get_w();
  for (size_t i = 0; i < U.N_c; i++) {
    u[i] *= Complex(z);
    v[i] *= Complex(z);
    w[i] *= Complex(z);
  }
  return _su3_accum(u, v, w);
}

inline _su3_accum operator*(const _su3_accum &U, const Complex &z) {
  return z * U;
}

template <> inline _su3_accum traceless_antiherm(const _su3_accum &x0) {
  _su3 Id; // by default is the identity
  _su3_accum x = x0;
  x = x - (trace(x) / double(x.N_c)) * Id;
  x = 0.5 * (x - x.dagger());
  return x;
}

using su3_accum = _su3_accum;

inline su3 accum_to_Group(const su3_accum &x) {
  su3 U(x.get_u(), x.get_v());
  U.restoreSU();
  return U;
}

/**
 * @brief accumulation type for SU(3)
 *
 * when summing SU(3) matrices the result is not unitary anymore. This type serves to
 * define the sum of the plaquettes appearing in the action and monomial forces
 * (covariant derivatives)
 *
 */
template <> struct accum_type<su3> {
  typedef su3_accum type;
};
