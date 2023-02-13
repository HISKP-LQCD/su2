// geometry.hh

#pragma once

#include <array>
#include <cstddef>
#include <iostream>
#include <numeric>

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

  inline void fatal_error(const std::string msg, char const *function_name) {
    std::cerr << "# FATAL ERROR.  " << function_name << "(): " << msg << "\nAborting.";
    std::cerr << "# Aborting\n";
    std::abort();
    return;
  }

} // namespace spacetime_lattice

/**
 * @brief standardized element access for a position-dependent quantity
 * 
 */
class geometry {
public:

  geometry(){}
  ~geometry(){}

  explicit geometry(const size_t _Lx,
                    const size_t _Ly,
                    const size_t _Lz,
                    const size_t _Lt)
    : Lx(_Lx), Ly(_Ly), Lz(_Lz), Lt(_Lt) {}
  size_t getLx() const { return (Lx); }
  size_t getLy() const { return (Ly); }
  size_t getLz() const { return (Lz); }
  size_t getLt() const { return (Lt); }
  size_t getIndex(const int t, const int x, const int y, const int z) const {
    size_t y0 = (t + Lt) % Lt;
    size_t y1 = (x + Lx) % Lx;
    size_t y2 = (y + Ly) % Ly;
    size_t y3 = (z + Lz) % Lz;
    return (((y0 * Lx + y1) * Ly + y2) * Lz + y3);
  }

private:
  size_t Lx, Ly, Lz, Lt;
};

/**
 * @brief standardized element access for a quantity like U_\mu(x)
 * 
 */
class geometry_link {
public:

  geometry_link(){}
  ~geometry_link(){}

  explicit geometry_link(const size_t _Lx,
                    const size_t _Ly,
                    const size_t _Lz,
                    const size_t _Lt,
                    const size_t _ndims)
    : Lx(_Lx), Ly(_Ly), Lz(_Lz), Lt(_Lt), ndims(_ndims) {}
  size_t getLx() const { return Lx; }
  size_t getLy() const { return Ly; }
  size_t getLz() const { return Lz; }
  size_t getLt() const { return Lt; }
  size_t getndims() const { return ndims; }
  size_t getIndex(const size_t t,
                  const size_t x,
                  const size_t y,
                  const size_t z,
                  const size_t _mu) const {
    size_t y0 = (t + Lt) % Lt;
    size_t y1 = (x + Lx) % Lx;
    size_t y2 = (y + Ly) % Ly;
    size_t y3 = (z + Lz) % Lz;
    size_t mu = (_mu + ndims) % ndims;
    return ((((y0 * Lx + y1) * Ly + y2) * Lz + y3) * ndims + mu);
  }

private:
  size_t Lx, Ly, Lz, Lt, ndims;
};

