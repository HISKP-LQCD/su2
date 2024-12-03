#pragma once

#include <random>

#include "gaugeconfig.hh"
#include "random_element.hh"
#include "su2.hh"

template <class T> void random_gauge_trafo(gaugeconfig<T> &U, const int seed) {
  std::mt19937 engine(seed);

  T rU, tmp;

  geometry Geom = U.get_geometry(); // geometry of the lattice
  const std::vector<size_t> L = Geom.get_L();
  const size_t n_dims = Geom.get_n_dims(); // number of dimensions
  const size_t N_pts = Geom.get_N_pts(); // number of dimensions
  for (size_t i = 0; i < N_pts; i++) {
    std::vector<size_t> x = spacetime_lattice::index_to_x<size_t>(i, L);
    std::vector<size_t> xminusmu = x;
    random_element(rU, engine, 1);
    for (size_t mu = 0; mu < n_dims; mu++) {
      U(x, mu) = rU * U(x, mu);

      xminusmu[mu] -= 1;
      U(xminusmu, mu) = U(xminusmu, mu) * rU.dagger();
      xminusmu[mu] += 1;
    }
  }

  return;
}
