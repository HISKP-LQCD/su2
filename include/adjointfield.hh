#pragma once

#include "geometry.hh"
#include <array>
#include <cassert>
#include <cmath>
#include <random>
#include <vector>

#include "adjoint_su2.hh"
#include "adjoint_su3.hh"
#include "adjoint_u1.hh"

// The following class will be used to deliver the
// adjoint type depending on the gauge group
template <typename Float, class Group> struct adjoint_type {
  typedef Group type;
};

template <typename Float> struct adjoint_type<Float, u1> {
  typedef adjointu1<Float> type;
};

template <typename Float> struct adjoint_type<Float, su2> {
  typedef adjointsu2<Float> type;
};

template <typename Float> struct adjoint_type<Float, su3> {
  typedef adjointsu3<Float> type;
};

/**
 * @brief field of adjoint matrices
 *
 * @tparam Float
 * @tparam Group
 */
template <typename Float, class Group> class adjointfield {
public:
  using value_type = typename adjoint_type<Float, Group>::type;
  template <class T> using nd_max_arr = typename spacetime_lattice::nd_max_arr<T>;

  adjointfield(const size_t Lx,
               const size_t Ly,
               const size_t Lz,
               const size_t Lt,
               const size_t ndims = 4)
    : Lx(Lx), Ly(Ly), Lz(Lz), Lt(Lt), volume(Lx * Ly * Lz * Lt), ndims(ndims) {
    data.resize(volume * 4);
  }
  adjointfield(const adjointfield &U)
    : Lx(U.getLx()),
      Ly(U.getLy()),
      Lz(U.getLz()),
      Lt(U.getLt()),
      volume(U.getVolume()),
      ndims(U.getndims()) {
    data.resize(volume * ndims);
    for (size_t i = 0; i < getSize(); i++) {
      data[i] = U[i];
    }
  }
  void flipsign() {
    for (size_t i = 0; i < getSize(); i++) {
      data[i].flipsign();
    }
  }
  size_t storage_size() const { return data.size() * sizeof(value_type); };
  size_t getLx() const { return (Lx); }
  size_t getLy() const { return (Ly); }
  size_t getLz() const { return (Lz); }
  size_t getLt() const { return (Lt); }
  size_t getndims() const { return (ndims); }
  size_t getVolume() const { return (volume); }
  size_t getSize() const { return (volume * ndims); }
  void operator=(const adjointfield<Float, Group> &U) {
    Lx = U.getLx();
    Ly = U.getLz();
    Lz = U.getLz();
    Lt = U.getLt();
    volume = U.getVolume();
    data.resize(U.getSize());
    for (size_t i = 0; i < U.getSize(); i++) {
      data[i] = U[i];
    }
  }

  value_type &operator()(
    size_t const t, size_t const x, size_t const y, size_t const z, size_t const mu) {
    return data[getIndex(t, x, y, z, mu)];
  }

  const value_type &operator()(size_t const t,
                               size_t const x,
                               size_t const y,
                               size_t const z,
                               size_t const mu) const {
    return data[getIndex(t, x, y, z, mu)];
  }

  value_type &operator()(std::vector<size_t> const &coords, size_t const mu) {
    return data[getIndex(coords[0], coords[1], coords[2], coords[3], mu)];
  }

  const value_type &operator()(std::vector<size_t> const &coords, size_t const mu) const {
    return data[getIndex(coords[0], coords[1], coords[2], coords[3], mu)];
  }

  template <class Type>
  value_type &operator()(const nd_max_arr<Type> &x, const size_t &mu) {
    return data[getIndex(x[0], x[1], x[2], x[3], mu)];
  }

  template <class Type>
  const value_type &operator()(const nd_max_arr<Type> &x, const size_t &mu) const {
    return data[getIndex(x[0], x[1], x[2], x[3], mu)];
  }

  value_type &operator[](size_t const index) { return data[index]; }

  const value_type &operator[](size_t const index) const { return data[index]; }

private:
  size_t Lx, Ly, Lz, Lt, volume, ndims;

  std::vector<value_type> data;

  size_t getIndex(const size_t t,
                  const size_t x,
                  const size_t y,
                  const size_t z,
                  const size_t mu) const {
    size_t y0 = (t + Lt) % Lt;
    size_t y1 = (x + Lx) % Lx;
    size_t y2 = (y + Ly) % Ly;
    size_t y3 = (z + Lz) % Lz;
    size_t _mu = (mu + ndims) % ndims;
    return ((((y0 * Lx + y1) * Ly + y2) * Lz + y3) * ndims + _mu);
  }
};

template <typename Float, class Group>
adjointfield<Float, su2> operator*(const Float &x, const adjointfield<Float, Group> &A) {
  adjointfield<Float, Group> res(A.getLx(), A.getLy(), A.getLz(), A.getLt());
  for (size_t i = 0; i < A.getSize(); i++) {
    res[i] = x * A[i];
  }
  return res;
}

template <typename Float>
Float operator*(const adjointfield<Float, _u1> &A, const adjointfield<Float, _u1> &B) {
  Float res = 0.;
  assert(A.getSize() == B.getSize());
  for (size_t i = 0; i < A.getSize(); i++) {
    res += A[i].geta() * B[i].geta();
  }
  return res;
}

template <typename Float>
Float operator*(const adjointfield<Float, su2> &A, const adjointfield<Float, su2> &B) {
  Float res = 0.;
  assert(A.getSize() == B.getSize());
  for (size_t i = 0; i < A.getSize(); i++) {
    res +=
      A[i].geta() * B[i].geta() + A[i].getb() * B[i].getb() + A[i].getc() * B[i].getc();
  }
  return res;
}

template <typename Float>
Float operator*(const adjointfield<Float, su3> &A, const adjointfield<Float, su3> &B) {
  Float res = 0.;
  assert(A.getSize() == B.getSize());
  for (size_t i = 0; i < A.getSize(); i++) {
    const std::array<Float, 8> arr_A = A[i].get_arr();
    const std::array<Float, 8> arr_B = B[i].get_arr();
    for (size_t k = 0; k < 8; k++) {
      res += arr_A[i] * arr_B[i];
    }
  }
  return res;
}

template <class URNG, typename Float>
void initnormal(URNG &engine, adjointfield<Float, _u1> &A) {
  std::normal_distribution<double> normal(0., 1.);
  for (size_t i = 0; i < A.getSize(); i++) {
    A[i].seta(Float(normal(engine)));
  }
  return;
}

template <class URNG, typename Float>
void initnormal(URNG &engine, adjointfield<Float, su2> &A) {
  std::normal_distribution<double> normal(0., 1.);
  for (size_t i = 0; i < A.getSize(); i++) {
    A[i].seta(Float(normal(engine)));
    A[i].setb(Float(normal(engine)));
    A[i].setc(Float(normal(engine)));
  }
  return;
}

template <class URNG, typename Float>
void initnormal(URNG &engine, adjointfield<Float, su3> &A) {
  std::normal_distribution<double> normal(0., 1.);
  for (size_t i = 0; i < A.getSize(); i++) {
    std::array<Float, 8> arr;
    for (size_t i = 0; i < 8; i++) {
      arr[i] = Float(normal(engine));
    }
    A[i].set_arr(arr);
  }
  return;
}

template <typename Float, class Group>
inline void zeroadjointfield(adjointfield<Float, Group> &A) {
  for (size_t i = 0; i < A.getSize(); i++) {
    A[i].setzero();
  }
}
