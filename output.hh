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

} // namespace output
