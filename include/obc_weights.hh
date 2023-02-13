/**
 * @file obc_weights.hh
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief weights of the action with open boundary conditions:
 * see eq. 2.7 of https://arxiv.org/pdf/1105.4749.pdf multipled by (-g_0^2) and
 * subtracted
 * @version 0.1
 * @date 2023-01-30
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "geometry.hh"

namespace obc {

  class weights {
  public:
    weights() {}
    ~weights() {}

    weights(const std::string &bc_type,
            const size_t Lx,
            const size_t Ly,
            const size_t Lz,
            const size_t Lt,
            const size_t ndims = spacetime_lattice::nd_max)
      : Lx(Lx), Ly(Ly), Lz(Lz), Lt(Lt), volume(Lx * Ly * Lz * Lt), ndims(ndims) {
      data.resize(volume, 1.0); // default weight is 1 (no obc)
      const geometry geom1(Lx, Ly, Lz, Lt);
      Geom = geom1;
      if (bc_type == "spatial_open") {
        this->apply_spatial_obc();
      } else {
        std::cerr << "Error, unsupported open boundary condition of type: " << bc_type << "\n";
        std::cerr << "Aborting\n";
        std::abort();
      }
    }

    std::vector<double> get_data() const { return data; }

    double &operator()(std::vector<size_t> const &coords) {
      return data[(*this).Geom.getIndex(coords[0], coords[1], coords[2], coords[3])];
    }

    double operator()(std::vector<size_t> const &coords) const {
      return data[(*this).Geom.getIndex(coords[0], coords[1], coords[2], coords[3])];
    }

  private:
    size_t Lx, Ly, Lz, Lt, volume, ndims;

    std::vector<double> data; // vector of values of the weights
    geometry Geom; // lattice geometry

    /**
     * @brief initialize the weights of eq. (2.7) of https://arxiv.org/pdf/1105.4749.pdf,
     * generalized to only-spatial open boundary conditions
     */
    void apply_spatial_obc() {
#pragma omp parallel for
      for (size_t x0 = 0; x0 < (*this).Lt; x0++) {
        for (size_t x1 = 0; x1 < (*this).Lx; x1++) {
          for (size_t x2 = 0; x2 < (*this).Ly; x2++) {
            for (size_t x3 = 0; x3 < (*this).Lz; x3++) {
              const std::vector<size_t> x = {x0, x1, x2, x3};

              const bool b1 = (x1 == 0 && ndims > 1);
              const bool b2 = (x2 == 0 && ndims > 2);
              const bool b3 = (x3 == 0 && ndims > 3);
              if (b1 || b2 || b3) {
                (*this)(x) = 0.0;
              }
            }
          }
        }
      }
    }
  };

} // namespace obc
