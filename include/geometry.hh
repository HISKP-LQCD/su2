// geometry.hh

#pragma once

#include <array>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <vector>

namespace spacetime_lattice {
  const size_t nd_max = 4; // maximum number of spacetime dimensions supported
  template <class T> using nd_max_arr = std::array<T, nd_max>;

  template <class T> size_t Npts_from_dims(const nd_max_arr<T> &dims) {
    size_t N = 1;
    for (size_t i = 0; i < nd_max; i++) {
      N *= dims[i];
    }
    return N;
  }

  /**
   * @brief number of positive-oriented L-shaped links starting from a point x
   * This function gives the following sum: $\sum_{\mu \neq \nu} 1.
   * The result is casted into double because it's often used in denominators
   * @param d number of dimensions
   * @return double value of the sum
   */
  inline double num_pLloops(const size_t &d) {
    return d * (d - 1);
  }

  /**
   * @brief num_pLloops(d)/2
   * This function gives the following sum: $\sum_{\mu < \nu} 1.
   * @param d
   * @return double
   */
  inline double num_pLloops_half(const int &d) {
    return num_pLloops(d) / 2.0;
  }

  /**
   * Transforming a vector x into a checkerboard index according to the equation (e.g.
   * n_dims=3): i = x0*L1*L2 + x1*L2 + x2
   *
   * NOTE: The routine is general, it works for all dimensions.
   *       The user needs to pass the vector of sized
   */
  inline size_t x_to_index(const std::vector<size_t> &x, const std::vector<size_t> &L) {
    size_t idx = x[0];
    size_t n_dims = L.size();
    for (size_t i = 1; i < n_dims; i++) {
      idx = L[i] * idx + x[i]; // at the end of the loop, idx has the correct expression
    }
    return idx;
  }

  /**
   * Transforming a checkerboard index to a vector "x",
   * where the convention is the same as for the inverse function x_to_index()
   *
   */
  inline std::vector<size_t> index_to_x(const size_t &i, std::vector<size_t> L) {
    size_t n_dims = L.size(); // number of dimensions
    std::vector<size_t> x(n_dims, 0);
    size_t i_sub = i;
    size_t den = 1;
    for (size_t mu = 0; mu < n_dims; mu++) {
      const size_t mu_rev = n_dims - mu - 1;
      x[mu_rev] = (i_sub % (den * L[mu_rev])) / den;
      i_sub -= den * x[mu_rev];
      den *= L[mu_rev];
    }
    return x;
  }

} // namespace spacetime_lattice

class geometry {
private:
  std::vector<size_t> L; // lattice sizes, e.g. Lt, Lx, Ly, Lz in 4 dimensions
  size_t n_dims = 0; // number of dimensions
  size_t N_pts = 0; // number of points of the lattice
  // std::vector<size_t> idx_x = {}; // checkerboard indices of "x"
  std::vector<std::vector<size_t>> idx_xplus_mu = {}; // checkerboard indices of x+\mu

  // set the arrays of indices for x + \mu, with periodic boundary conditions
  void set_xpmu_idx_pbc() {
    idx_xplus_mu.resize(n_dims);
    for (size_t i = 0; i < N_pts; i++) {
      std::vector<size_t> x = spacetime_lattice::index_to_x(i, (*this).L);
      // idx_x.push_back(i);
      for (size_t mu = 0; mu < n_dims; mu++) {
        x[mu] = (x[mu] + 1) % L[mu];
        idx_xplus_mu[mu].push_back(spacetime_lattice::x_to_index(x, L));
        x[mu] = (x[mu] - 1 + L[mu]) % L[mu];
      }
    }
  }

public:
  geometry() {}
  ~geometry() {}

  explicit geometry(const size_t _Lx,
                    const size_t _Ly,
                    const size_t _Lz,
                    const size_t _Lt) {
    L = {_Lt, _Lx, _Lz, _Lt};
    n_dims = L.size();
    idx_xplus_mu.resize(n_dims);
    N_pts = std::accumulate(L.begin(), L.end(), 1.0, std::multiplies<double>());
  }

  explicit geometry(const std::vector<size_t> &_L) {
    L = _L;
    n_dims = L.size();
    N_pts = std::accumulate(L.begin(), L.end(), 1.0, std::multiplies<double>());
    this->set_xpmu_idx_pbc();
  }

  size_t getLt() const { return L[0]; }
  size_t getLx() const { return L[1]; }
  size_t getLy() const { return L[2]; }
  size_t getLz() const { return L[3]; }

  size_t get_n_dims() const { return L.size(); }

  size_t get_N_pts() const { return N_pts; }
  // std::vector<size_t> get_idx_x() const { return idx_x; }
  std::vector<std::vector<size_t>> get_idx_xplus_mu() const { return idx_xplus_mu; }

  size_t getIndex(const int t, const int x, const int y, const int z) const {
    size_t Lt = L[0];
    size_t Lx = L[1];
    size_t Ly = L[2];
    size_t Lz = L[3];
    size_t y0 = (t + Lt) % Lt;
    size_t y1 = (x + Lx) % Lx;
    size_t y2 = (y + Ly) % Ly;
    size_t y3 = (z + Lz) % Lz;
    return (((y0 * Lx + y1) * Ly + y2) * Lz + y3);
  }

};
