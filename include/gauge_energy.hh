/**
 * @file gauge_energy.hpp
 * @author simone-romiti (simone.romiti@uni-bonn.de)
 * @brief
 * @version 0.1
 * @date 2022-05-20
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#include "metric.hpp"

#include "flat-gauge_energy.hpp"
#include "rotating-gauge_energy.hpp"

template <class T, class M>
double gauge_energy_g(const M& g_munu, const gaugeconfig<T> &U, bool spatial_only = false);

template <class T>
double gauge_energy_g<T, metric::flat>(const gaugeconfig<T> &U, bool spatial_only = false) {
  return flat_spacetime::gauge_energy(U, spatial_only);
}

template <class T>
double gauge_energy(const gaugeconfig<T> &U, bool spatial_only = false) {
  return gauge_energy_g(U, spatial_only);
}
