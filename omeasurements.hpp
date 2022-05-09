// omeasurements.hpp
/**
 * @brief routines for online/offline measurements of various observables.
 * This file defines functions which handle the measure and printing of observables over a
 * given gauge configuration.
 */

#pragma once

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "wilsonloop.hh"
#include"gradient_flow.hh"

namespace omeasurements {

  /**
   * @brief compute and print the wilson loop of a given configuration
   *
   * @tparam Group
   * @param U gauge config
   * @param i configurationn index
   * @param confdir
   */
  template <class Group>
  void
  meas_wilson_loop(const gaugeconfig<Group> &U, const size_t &i, const std::string &confdir) {
    std::ostringstream os;
    os << confdir + "/wilsonloop.";
    auto prevw = os.width(6);
    auto prevf = os.fill('0');
    os << i;
    os.width(prevw);
    os.fill(prevf);
    os << ".dat" << std::ends;
    compute_all_loops(U, os.str());

    return;
  }

  /**
   * @brief compute and print the gradient flow of a given configuration
   *
   * @tparam Group
   * @param U gauge config
   * @param i configurationn index
   * @param confdir
   */
  template <class Group>
  void meas_gradient_flow(const gaugeconfig<Group> &U,
                     const size_t &i,
                     const std::string &confdir,
                     const double &tmax) {
    std::ostringstream os;
    os << confdir + "/gradient_flow.";
    auto prevw = os.width(6);
    auto prevf = os.fill('0');
    os << i;
    os.width(prevw);
    os.fill(prevf);
    gradient_flow(U, os.str(), tmax);

    return;
  }

} // namespace omeasurements