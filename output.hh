// output.hh
/**
 * @file output.hh
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief routines for outputting data
 * @version 0.1
 * @date 2022-05-09
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <iomanip>
#include <sstream>
#include <string>

#include "parameters.hh"

namespace output {
  namespace gp = global_parameters;

  /**
   * @brief Get the beginning of the configuration path's string
   * This function returns the sub-string containing the path to the configuration,
   * without the information about the index of the sweep(Metropolis)/trajectory(HMC),etc.
   * It is assumed that 'sparams' has the attributes:
   * - conf_dir
   * - conf_basename
   * - beta_str_width
   * @tparam S type of structure containing the conf_basename attribute
   * @param pparams physical parameters
   * @param sparams parameters specific to the given
   * @return std::string
   */
  template <class S>
  std::string get_conf_path_basename(const gp::physics &pparams, const S &sparams) {
    // set up name for configs
    const std::string conf_basename = sparams.conf_basename;
    std::stringstream ss;
    ss << sparams.conf_basename << ".";
    ss << pparams.Lx << "." << pparams.Ly << "." << pparams.Lz << "." << pparams.Lt;
    if (pparams.rotating_frame) {
      ss << ".Omega_" << pparams.Omega;
    }
    ss << ".b" << std::fixed << std::setprecision(sparams.beta_str_width) << pparams.beta;
    if (pparams.anisotropic) {
      ss << ".x" << std::fixed << std::setprecision(sparams.beta_str_width) << pparams.xi;
    }

    return sparams.conf_dir+ss.str();
  }

  std::string get_filename_fine(const gp::physics &pparams,
                                const gp::measure_u1 &mparams) {
    std::ostringstream filename_fine;

    filename_fine << mparams.resdir << "/"
                  << "result" << pparams.ndims - 1 << "p1d.u1potential.rotated.Nt"
                  << pparams.Lt << ".Ns" << pparams.Lx << ".b" << std::fixed
                  << std::setprecision(mparams.beta_str_width) << pparams.beta << ".xi"
                  << std::fixed << std::setprecision(mparams.beta_str_width) << pparams.xi
                  << ".nape" << mparams.n_apesmear << ".alpha" << std::fixed
                  << mparams.alpha << "finedistance" << std::ends;

    return filename_fine.str();
  }

  std::string get_filename_coarse(const gp::physics &pparams,
                                  const gp::measure_u1 &mparams) {
    std::ostringstream f;

    f << mparams.resdir << "/"
      << "result" << pparams.ndims - 1 << "p1d.u1potential.rotated.Nt" << pparams.Lt
      << ".Ns" << pparams.Lx << ".b" << std::fixed
      << std::setprecision(mparams.beta_str_width) << pparams.beta << ".xi" << std::fixed
      << std::setprecision(mparams.beta_str_width) << pparams.xi << ".nape"
      << mparams.n_apesmear << ".alpha" << std::fixed << mparams.alpha << "coarsedistance"
      << std::ends;

    return f.str();
  }

  std::string get_filename_nonplanar(const gp::physics &pparams,
                                     const gp::measure_u1 &mparams) {
    std::ostringstream f;

    f << mparams.resdir << "/"
      << "result" << pparams.ndims - 1 << "p1d.u1potential.Nt" << pparams.Lt << ".Ns"
      << pparams.Lx << ".b" << std::fixed << std::setprecision(mparams.beta_str_width)
      << pparams.beta << ".xi" << std::fixed << std::setprecision(mparams.beta_str_width)
      << pparams.xi << ".nape" << mparams.n_apesmear << ".alpha" << std::fixed
      << mparams.alpha << "nonplanar" << std::ends;

    return f.str();
  }

  namespace hmc {

    std::string get_header(const std::string &sep = " ") {
      std::stringstream ss; // header: column names in the output
      ss << "i" << sep << "getaccept" << sep << "E*A" << sep << "dH" << sep << "rho"
         << sep << "ddH" << sep << "Q"
         << "\n";
      return ss.str();
    }

  } // namespace hmc

} // namespace output
