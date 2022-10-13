// io.hh
/**
 * @file io.hh
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief routines for input-output of data
 * @version 0.1
 * @date 2022-05-09
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include "parameters.hh"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <xtensor/xcsv.hpp>

namespace io {
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
    ss << sparams.conf_basename;
    if (sparams.lenghty_conf_name) {
      ss << "." << pparams.Lx << "." << pparams.Ly << "." << pparams.Lz << "."
         << pparams.Lt;
      if (pparams.rotating_frame) {
        ss << ".Omega_" << pparams.Omega;
      }
      ss << ".b" << std::fixed << std::setprecision(sparams.beta_str_width)
         << pparams.beta;
      if (pparams.anisotropic) {
        ss << ".x" << std::fixed << std::setprecision(sparams.beta_str_width)
           << pparams.xi;
      }
    }

    return sparams.conf_dir + "/" + ss.str();
  }

  /**
   * @brief print the content of an xtensor expression to file
   * This functions is needed because the current (October 2022) master branch of xtensor
   * doesn't support a custom delimiter (and neither a header) for dumping to file
   * @tparam E expression type
   */
  template <class E>
  void xtensor_to_stream(std::ofstream &stream,
                         const xt::xexpression<E> &e,
                         const std::string &sep = " ",
                         const std::string &h = "") {
    using size_type = typename E::size_type;
    const E &ex = e.derived_cast();
    if (ex.dimension() != 2) {
      XTENSOR_THROW(std::runtime_error, "Only 2-D expressions can be serialized to CSV");
    }
    if (h != "") {
      stream << h << std::endl;
    }
    size_type nbrows = ex.shape()[0], nbcols = ex.shape()[1];
    auto st = ex.stepper_begin(ex.shape());
    for (size_type r = 0; r != nbrows; ++r) {
      for (size_type c = 0; c != nbcols; ++c) {
        stream << *st;
        if (c != nbcols - 1) {
          st.step(1);
          stream << sep;
        } else {
          st.reset(1);
          st.step(0);
          stream << std::endl;
        }
      }
    }
    return;
  }

  namespace measure {

    /**
     * gets filename with several physical characteristics that are contained in pparams
     * and S, which should be either metropolis_u1 or measure_u1
     * **/
    template <class S>
    std::string get_filename_fine(const gp::physics &pparams, const S &mparams) {
      std::ostringstream filename_fine;

      filename_fine << mparams.res_dir << "/"
                    << "result" << pparams.ndims - 1 << "p1d.u1potential.rotated.Nt"
                    << pparams.Lt << ".Ns" << pparams.Lx << ".b" << std::fixed
                    << std::setprecision(mparams.beta_str_width) << pparams.beta << ".xi"
                    << std::fixed << std::setprecision(mparams.beta_str_width)
                    << pparams.xi << ".nape" << mparams.n_apesmear << ".alpha"
                    << std::fixed << mparams.alpha << "finedistance" << std::ends;

      return filename_fine.str();
    }

    /**
     * gets filename with several physical characteristics that are contained in pparams
     * and S, which should be either metropolis_u1 or measure_u1
     * **/
    template <class S>
    std::string get_filename_coarse(const gp::physics &pparams, const S &mparams) {
      std::ostringstream f;

      f << mparams.res_dir << "/"
        << "result" << pparams.ndims - 1 << "p1d.u1potential.rotated.Nt" << pparams.Lt
        << ".Ns" << pparams.Lx << ".b" << std::fixed
        << std::setprecision(mparams.beta_str_width) << pparams.beta << ".xi"
        << std::fixed << std::setprecision(mparams.beta_str_width) << pparams.xi
        << ".nape" << mparams.n_apesmear << ".alpha" << std::fixed << mparams.alpha
        << "coarsedistance" << std::ends;

      return f.str();
    }

    /**
     * gets filename with several physical characteristics that are contained in pparams
     * and S, which should be either metropolis_u1 or measure_u1
     * **/
    template <class S>
    std::string get_filename_nonplanar(const gp::physics &pparams, const S &mparams) {
      std::ostringstream f;

      f << mparams.res_dir << "/"
        << "result" << pparams.ndims - 1 << "p1d.u1potential.Nt" << pparams.Lt << ".Ns"
        << pparams.Lx << ".b" << std::fixed << std::setprecision(mparams.beta_str_width)
        << pparams.beta << ".xi" << std::fixed
        << std::setprecision(mparams.beta_str_width) << pparams.xi << ".nape"
        << mparams.n_apesmear << ".alpha" << std::fixed << mparams.alpha << "nonplanar"
        << std::ends;

      return f.str();
    }

    /**
     * @brief writes the headers for the results of the planar Wilson-Loops
     * in the files. coarse is the spacial direction, fine the temporal direction.
     * a_t<=a_s.
     * @param pparams physical parameters
     * @param mparams parameters specific to measurement. May change if a non-valid
     * dimension is given, so cannot be const. S is metropolis_u1 or measure_u1.
     * @param filenamecoarse, filenamefine: names of the files where results are stored.
     * **/
    template <class S>
    void set_header_planar(const gp::physics &pparams,
                           S &mparams,
                           const std::string &filename_coarse,
                           const std::string &filename_fine) {
      //~ open file for saving results
      std::ofstream resultfile;
      if (pparams.ndims == 2) {
        std::cerr << "Currently not working for dim = 2, no measurements for the "
                     "potential will be made"
                  << std::endl;
        mparams.potentialplanar = false;
      }

      //~ print heads of columns: W(r, t), W(x, y)
      if (!mparams.append && (pparams.ndims == 3 || pparams.ndims == 4)) {
        resultfile.open(filename_fine, std::ios::out);
        resultfile << "## ";
        for (size_t t = 1; t <= pparams.Lt * mparams.sizeWloops; t++) {
          for (size_t x = 1; x <= pparams.Lx * mparams.sizeWloops; x++) {
            resultfile << "W(x=" << x << ",t=" << t << ",y=" << 0 << ")  ";
          }
        }
        resultfile << "counter";
        resultfile << std::endl;
        resultfile.close();

        resultfile.open(filename_coarse, std::ios::out);
        resultfile << "## ";
        for (size_t y = 1; y <= pparams.Ly * mparams.sizeWloops; y++) {
          for (size_t x = 1; x <= pparams.Lx * mparams.sizeWloops; x++) {
            resultfile << "W(x=" << x << ",t=" << 0 << ",y=" << y << ")  ";
          }
        }
        resultfile << "counter";
        resultfile << std::endl;
        resultfile.close();
      }
    }

    /**
     * @brief writes the headers for the results of the non-planar Wilson-Loops in the
     * file.
     * @param pparams physical parameters
     * @param mparams parameters specific to measurement. May change if a non-valid
     * dimension is given, so cannot be const. S is metropolis_u1 or measure_u1.
     * @param filenamenonplanar: name of the files where results are stored.
     * **/
    template <class S>
    void set_header_nonplanar(const gp::physics &pparams,
                              S &mparams,
                              const std::string &filename_nonplanar) {
      //~ open file for saving results
      std::ofstream resultfile;
      size_t maxsizenonplanar = (pparams.Lx < 4) ? pparams.Lx : 4;

      if (pparams.ndims == 2 || pparams.ndims == 4) {
        std::cerr << "Currently not working for dim = 2 and dim = 4, no nonplanar "
                     "measurements will be made"
                  << std::endl;
        mparams.potentialnonplanar = false;
      }

      //~ print heads of columns
      if (!mparams.append && (pparams.ndims == 3)) {
        resultfile.open(filename_nonplanar, std::ios::out);
        resultfile << "## ";
        for (size_t t = 0; t <= pparams.Lt * mparams.sizeWloops; t++) {
          for (size_t x = 0; x <= maxsizenonplanar; x++) {
            for (size_t y = 0; y <= maxsizenonplanar; y++) {
              resultfile << "W(x=" << x << ",t=" << t << ",y=" << y << ")  ";
            }
          }
        }
        resultfile << "counter";
        resultfile << std::endl;
        resultfile.close();
      }
    }
  } // namespace measure

  namespace hmc {

    std::string g_nconf_counter = "nconf_counter.txt"; // global variable: name of file

    std::string get_header(const std::string &sep = " ") {
      std::stringstream ss; // header: column names in the io
      ss << "i" << sep << "getaccept" << sep << "E*A" << sep << "dH" << sep << "rho"
         << sep << "ddH" << sep << "Q"
         << "\n";
      return ss.str();
    }

    /**
     * @brief index and path of the last saved configuration
     *
     * @param conf_dir directory with nconf_counter.txt inside
     * @return std::array<std::string, 2> {"i","/path/to/conf.i"}
     */
    std::vector<std::string> read_nconf_counter(const std::string &conf_dir) {
      const std::string input_file = conf_dir + g_nconf_counter;

      if (!boost::filesystem::exists(input_file)) {
        std::cerr << "Error from " << __func__ << "\n";
        std::cerr << input_file << ": no such file or directory. Aborting.\n";
        std::abort();
      }

      std::ifstream nconf_counter(input_file);
      std::string line;
      std::getline(nconf_counter, line); // 1st line is just a header
      std::getline(nconf_counter, line);
      nconf_counter.close();

      std::vector<std::string> v;
      boost::split(v, line, boost::is_any_of(" "));
      return {v[0], v[1], v[2]};
    }

    /**
     * @brief saving the index and path of the last gauge configuration
     *
     * @param conf_dir directory containing the configurations
     * @param i trajectory index of the last configuration
     * @param path_conf full path to the last configuration
     */
    void update_nconf_counter(const std::string &conf_dir,
                              const double &heat,
                              const size_t &i,
                              const std::string &path_conf) {
      // saving index of the last configuration
      std::ofstream nconf_counter;
      nconf_counter.open(conf_dir + g_nconf_counter, std::ios::out);
      nconf_counter << "heat i path_conf\n";
      nconf_counter << heat << " " << i << " " << path_conf;
      nconf_counter.close();
    }

  } // namespace hmc

} // namespace io
