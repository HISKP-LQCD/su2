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

#include "flat-gradient_flow.hh"
#include "links.hpp"
#include "operators.hpp"
#include "parameters.hh"
#include "propagator.hpp"
#include "smearape.hh"
#include "wilsonloop.hh"
#include <boost/filesystem.hpp>

namespace omeasurements {

  namespace fsys = boost::filesystem;

  /**
   * @brief compute and print the wilson loop of a given configuration
   *
   * @tparam Group
   * @param U gauge config
   * @param i configurationn index
   * @param conf_dir
   */
  template <class Group>
  void meas_wilson_loop(const gaugeconfig<Group> &U,
                        const size_t &i,
                        const std::string &conf_dir) {
    std::ostringstream os;
    os << conf_dir + "/wilsonloop.";
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
   * @tparam sparams struct containing info on computation and output
   * @param U gauge config
   * @param i configurationn index
   * @param conf_dir
   */
  template <class Group, class sparams>
  void
  meas_gradient_flow(const gaugeconfig<Group> &U, const size_t &i, const sparams &S) {
    const std::string res_dir = S.res_dir;
    const double eps = S.epsilon_gradient_flow;
    const double tmax = S.tmax;

    std::ostringstream os;
    os << res_dir + "/gradient_flow.";
    auto prevw = os.width(6);
    auto prevf = os.fill('0');
    os << i;
    os.width(prevw);
    os.fill(prevf);
    flat_spacetime::gradient_flow(U, os.str(), tmax, eps);

    return;
  }

  /**
   * @brief measure and print the (staggered) pion correlator
   * It is assumed that sparams (specific parameters) contain the following attributes:
   * -
   * @tparam Group
   * @tparam sparams struct containing info on computation and output
   * @param U gauge configuration
   * @param i trajectory index
   * @param S
   */
  template <class Group, class sparams>
  void meas_pion_correlator(const gaugeconfig<Group> &U,
                            const size_t &i,
                            const double &m,
                            const sparams &S) {
    std::ostringstream oss;
    oss << S.res_dir + "/C_pion.";
    auto prevw = oss.width(6);
    auto prevf = oss.fill('0');
    oss << i;
    oss.width(prevw);
    oss.fill(prevf);

    const std::string path = oss.str();
    std::ofstream ofs(path, std::ios::out);

    ofs << "t C(t)\n";
    const std::vector<double> Cpi = staggered::C_pion<Group>(
      U, m, S.solver, S.tolerance_cg, S.solver_verbosity, S.seed_pf);
    for (size_t i = 0; i < U.getLt(); i++) {
      ofs << std::scientific << std::setprecision(16) << i << " " << Cpi[i] << "\n";
    }
    ofs.close();

    return;
  }

  /**
   * @brief measure the glueball 0^{PC} correlators 
   * the interpolators are average plaquettes: sum_{\nu<\mu} U_{\mu \nu} 
   * measure the glueball correlators 0^{++}, 0^{+-}, 0^{-+}, 0^{--}
   * @tparam Group
   * @tparam sparams struct containing info on computation and output
   * @param U gauge configuration
   * @param i trajectory index
   * @param S specific parameters
   */
  template <class Group, class sparams>
  void meas_glueball_correlator_U_munu(const gaugeconfig<Group> &U0,
                                             const size_t &i,
                                             const sparams &S, const  bool& spatial_only) {
    typedef typename accum_type<Group>::type accum;

    std::ostringstream oss;
    std::string subdir = "U_munu/";
    if(spatial_only){subdir="U_ij/";}
    oss << S.res_dir + "/"+subdir+"C_glueball";
    fsys::create_directories(fsys::absolute(S.res_dir + "/" + subdir + "/"));

    gaugeconfig<Group> U = U0;
    if (S.glueball.doAPEsmear) {
      for (size_t i = 0; i < S.glueball.nAPEsmear; i++) {
        smearlatticeape<Group>(U, S.glueball.alphaAPEsmear, true);
      }
      oss << "_smearAPEn" << S.glueball.nAPEsmear << "alpha"
          << S.glueball.alphaAPEsmear;
    }
    oss << "_";
    auto prevw = oss.width(8);
    auto prevf = oss.fill('0');
    oss << i;
    oss.width(prevw);
    oss.fill(prevf);

    const std::string path = oss.str();
    std::ofstream ofs(path, std::ios::out);

    std::vector<double> sinks[4]; // sink at 't'
    // std::vector<double> sources[4]; // source at 't'
    for (size_t i_PC = 0; i_PC < 4; i_PC++) {
      sinks[i_PC].resize(U.getLt());
      // sources[i_PC].resize(U.getLt());
    }

    for (size_t t = 0; t < U.getLt(); t++) {
      const accum Uij = operators::get_rest_tr_sum_U_munu<accum, Group>(U, t, spatial_only, false);
      const accum PUij = operators::get_rest_tr_sum_U_munu<accum, Group>(U, t, spatial_only, true);

      size_t i_PC = 0;
      for (int sP = 1; sP >= -1; sP -= 2) {
        for (int iC = 1; iC >= 0; iC--) {
          const bool C = bool(iC);
          std::complex<double> comb_snk = Uij + double(sP) * PUij;
          comb_snk /= 2.0;

          double snk;
          if (C) {
            snk = std::real(comb_snk);
          } else {
            snk = std::imag(comb_snk);
          }
          // snk = 1.0 - snk; // vacuum subtraction
          sinks[i_PC][t] = snk;

          // std::complex<double> comb_src =
          //   operators::get_tr_sum_U_ij<accum, Group>(U, {int(t), 0, 0, 0}, false) +
          //   double(sP) *
          //     operators::get_tr_sum_U_ij<accum, Group>(U, {int(t), 0, 0, 0}, true);
          // comb_src /= 2.0;

          // double src;
          // if (C) {
          //   src = std::real(comb_src);
          // } else {
          //   src = std::imag(comb_src);
          // }

          // sources[i_PC][t] = 1.0 - src; // vacuum subtraction

          ++i_PC;
        }
      }
    }

    ofs << "t C_{++}(t) C_{+-}(t) C_{-+}(t) C_{--}(t)" << std::endl; // header
    const size_t T_ext = U.getLt(); // lattice temporal time extent
    for (size_t t = 0; t < T_ext; t++) {
      ofs << t;
      for (size_t i_PC = 0; i_PC < 4; i_PC++) {
        double Ct = 0.0;
        for (size_t tau = 0; tau < T_ext; tau++) {
          Ct += sinks[i_PC][(t + tau) % T_ext] * sinks[i_PC][tau];
        }
        Ct /= double(T_ext); // average over all times
        ofs << " " << std::scientific << std::setprecision(16) << Ct;
      }
      ofs << std::endl;
    }

    ofs.close();

    return;
  }

  /**
   * @brief measure the glueball 0^{PC} correlators for the GEVP
   * measure the glueball correlators ij for the 0^{++}, 0^{+-}, 0^{-+}, 0^{--} glueballs.
   * j<i and i runs from 0 to L-1, where L is the maximum lenght of the loops specified in
   * the 'sparams' struct.
   * @tparam Group
   * @tparam sparams struct containing info on computation and output
   * @param U gauge configuration
   * @param i trajectory index
   * @param S specific parameters
   */
  template <class Group, class sparams>
  void meas_glueball_correlator_GEVP(const gaugeconfig<Group> &U0,
                                     const size_t &i,
                                     const sparams &S) {
    typedef typename accum_type<Group>::type accum;

    std::ostringstream oss_dir, oss_name;
    oss_dir << S.res_dir + "/";

    std::ostringstream oss_details;
    oss_details << "smearAPEn" << S.glueball.nAPEsmear << "alpha"
                << S.glueball.alphaAPEsmear;
    std::string meas_details = oss_details.str();

    if (S.glueball.use_res_dir) {
      oss_dir << meas_details + "/";
    }

    // oss_name << "C_glueball";

    gaugeconfig<Group> U = U0;
    if (S.glueball.doAPEsmear) {
      for (size_t i = 0; i < S.glueball.nAPEsmear; i++) {
        smearlatticeape<Group>(U, S.glueball.alphaAPEsmear, true);
      }
      if (S.glueball.lengthy_file_name) {
        oss_name << "_" << meas_details;
      }
    }
    oss_name << "_";
    auto prevw = oss_name.width(8);
    auto prevf = oss_name.fill('0');
    oss_name << i;
    oss_name.width(prevw);
    oss_name.fill(prevf);

    const size_t max_length_loops = S.glueball.max_length_loops;
    std::vector<links::path> clpaths(0);

    for (size_t ll = 1; ll <= max_length_loops; ll++) {
      const links::closed_paths clpf(U.getndims(), 0, ll);
      const std::vector<links::path> cll =
        clpf.get_paths(); // vector of closed paths found
      clpaths.insert(clpaths.end(), cll.begin(), cll.end()); // all closed paths
    }

    const size_t nl = clpaths.size();

    // phi_i(t)^{PC}
    std::vector<std::vector<std::array<std::array<double, 2>, 2>>> phi(nl);

    for (size_t i = 0; i < nl; i++) {
      phi[i].resize(U.getLt());
      for (size_t t = 0; t < U.getLt(); t++) {
        const std::complex<double> Pp =
          operators::trace_rest_loop<std::complex<double>, Group>(t, U, clpaths[i],
                                                                  false);
        const std::complex<double> Pm =
          operators::trace_rest_loop<std::complex<double>, Group>(t, U, clpaths[i], true);
        phi[i][t][0][0] = (Pp + Pm).real(); // PC=++
        phi[i][t][0][1] = (Pp + Pm).imag(); // PC=+-
        phi[i][t][1][0] = (Pp - Pm).real(); // PC=-+
        phi[i][t][1][1] = (Pp - Pm).imag(); // PC=--
      }
    }

    const size_t T_ext = U.getLt(); // lattice temporal time extent
    for (size_t i = 0; i < nl; i++) {
      for (size_t j = 0; j <= i; j++) { // C_{ij} == C_{ji}

        const std::string dir_ij = oss_dir.str() + std::to_string(i) + "_" +
                                   std::to_string(j) + "/"; // directory path
        fsys::create_directories(fsys::absolute(dir_ij)); // creating directory
        const std::string path =
          dir_ij + "C_glueball" + oss_name.str(); // full path of output file

        std::ostringstream oss_ij;
        oss_ij << "t C_{++}(t) C_{+-}(t) C_{-+}(t) C_{--}(t)" << std::endl; // header
        for (size_t t = 0; t < T_ext; t++) {
          oss_ij << t;
          for (size_t P = 0; P <= 1; P++) {
            for (size_t C = 0; C <= 1; C++) {
              double Ct = 0.0;
              for (size_t tau = 0; tau < T_ext; tau++) {
                Ct += phi[i][(t + tau) % T_ext][P][C] * phi[j][tau][P][C];
              }
              Ct /= double(T_ext); // average over all times
              oss_ij << " " << std::scientific << std::setprecision(16) << Ct;
            }
          }
          oss_ij << std::endl;
        }

        std::ofstream ofs(path, std::ios::out);
        ofs << oss_ij.str();
        ofs.close();
      }
    }

    return;
  }

  /**
   * measures the planar wilson loops in temporal and spacial direction
   * -(0,1) loops for temporal, (1,2) and (1,3) loops for spacial
   * writes one line per configuration into the resultfiles.
   * At the moment, measurements are implemented for dim=3,4.
   * @param U holds the gauge-configuration whose loops are measured
   * @param pparams holds information about the size of the lattice
   * @param sizeWloops: maximum extent up to which loops are measured
   * @param filenames: files into which the results of the measurements are written
   * @param i: index of the configuration
   * **/
  template <class Group>
  void meas_loops_planar_pot(const gaugeconfig<Group> &U,
                             const global_parameters::physics &pparams,
                             const size_t &sizeWloops,
                             const std::string &filename_coarse,
                             const std::string &filename_fine,
                             const size_t &i) {
    double loop;
    std::ofstream resultfile;
    //~ //calculate wilsonloops for potential
    if (pparams.ndims == 4) {
      resultfile.open(filename_fine, std::ios::app);
      for (size_t t = 1; t <= pparams.Lt * sizeWloops; t++) {
        for (size_t x = 1; x <= pparams.Lx * sizeWloops; x++) {
          loop = wilsonloop_non_planar(U, {t, x, 0, 0});
          resultfile << std::setw(14) << std::scientific << loop / U.getVolume() << "  ";
        }
      }
      resultfile << i;
      resultfile << std::endl;
      resultfile.close();

      resultfile.open(filename_coarse, std::ios::app);
      for (size_t y = 1; y <= pparams.Ly * sizeWloops; y++) {
        for (size_t x = 1; x <= pparams.Lx * sizeWloops; x++) {
          loop = wilsonloop_non_planar(U, {0, x, y, 0});
          loop += wilsonloop_non_planar(U, {0, x, 0, y});
          resultfile << std::setw(14) << std::scientific << loop / U.getVolume() / 2.0
                     << "  ";
        }
      }
      resultfile << i;
      resultfile << std::endl;
      resultfile.close();
    }
    if (pparams.ndims == 3) {
      resultfile.open(filename_fine, std::ios::app);
      for (size_t t = 1; t <= pparams.Lt * sizeWloops; t++) {
        for (size_t x = 1; x <= pparams.Lx * sizeWloops; x++) {
          loop = wilsonloop_non_planar(U, {t, x, 0});
          //~ loop  += wilsonloop_non_planar(U, {t, 0, x});
          resultfile << std::setw(14) << std::scientific << loop / U.getVolume() << "  ";
        }
      }
      resultfile << i;
      resultfile << std::endl;
      resultfile.close();

      resultfile.open(filename_coarse, std::ios::app);
      for (size_t y = 1; y <= pparams.Ly * sizeWloops; y++) {
        for (size_t x = 1; x <= pparams.Lx * sizeWloops; x++) {
          loop = wilsonloop_non_planar(U, {0, x, y});
          //~ loop += wilsonloop_non_planar(U, {0, y, x});
          resultfile << std::setw(14) << std::scientific << loop / U.getVolume() << "  ";
        }
      }
      resultfile << i;
      resultfile << std::endl;
      resultfile.close();
    }
  }

  /**
   * measures the nonplanar wilson loops in temporal and spacial direction
   * writes one line per configuration into the resultfiles.
   * At the moment, measurements are implemented for dim=3.
   * all possible loops (t,x,y) with x,y <=min(4, Lx), t<Lt*sizeWloops are measured.
   * @param U holds the gauge-configuration whose loops are measured
   * @param pparams holds information about the size of the lattice
   * @param sizeWloops: maximum extent up to which loops are measured
   * @param filenames: files into which the results of the measurements are written
   * @param i: index of the configuration
   * **/
  template <class Group>
  void meas_loops_nonplanar_pot(const gaugeconfig<Group> &U,
                                const global_parameters::physics &pparams,
                                const size_t &sizeWloops,
                                const std::string &filename_nonplanar,
                                const size_t &i) {
    double loop;
    std::ofstream resultfile;
    size_t maxsizenonplanar = (pparams.Lx < 4) ? pparams.Lx : 4;
    resultfile.open(filename_nonplanar, std::ios::app);
    for (size_t t = 0; t <= pparams.Lt * sizeWloops; t++) {
      for (size_t x = 0; x <= maxsizenonplanar; x++) {
        for (size_t y = 0; y <= maxsizenonplanar; y++) {
          loop = wilsonloop_non_planar(U, {t, x, y});
          resultfile << std::setw(14) << std::scientific << loop / U.getVolume() << "  ";
        }
      }
    }
    resultfile << i;
    resultfile << std::endl;
    resultfile.close();
  }

} // namespace omeasurements
