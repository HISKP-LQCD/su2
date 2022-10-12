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
#include "glueballs.hpp"
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
                        const std::string &res_dir) {
    std::ostringstream os;
    os << res_dir + "/wilsonloop.";
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
  void meas_gradient_flow(const gaugeconfig<Group> &U,
                          const size_t &i,
                          const global_parameters::physics &pparams,
                          const sparams &S) {
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
    flat_spacetime::gradient_flow(U, os.str(), tmax, eps, pparams.xi);

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

  namespace from_smeared_field {

    // /**
    //  * @brief measure the glueball 0^{PC} correlators
    //  * the interpolators are average plaquettes: sum_{\nu<\mu} U_{\mu \nu}
    //  * measure the glueball correlators 0^{++}, 0^{+-}, 0^{-+}, 0^{--}
    //  * @tparam Group
    //  * @tparam sparams struct containing info on computation and output
    //  * @param U gauge configuration
    //  * @param i trajectory index
    //  * @param S specific parameters
    //  */
    // template <class Group, class sparams>
    // void meas_glueball_correlator_U_munu(const gaugeconfig<Group> &U,
    //                                      const size_t &i,
    //                                      const size_t &nAPEsmear,
    //                                      const sparams &S,
    //                                      const bool &spatial_only) {
    //   typedef typename accum_type<Group>::type accum;

    //   std::ostringstream oss;
    //   oss << S.res_dir + "/";

    //   // apply smearing
    //   if (S.glueball.doAPEsmear) {
    //     oss << "/smearAPEn" << nAPEsmear << "alpha" << S.glueball.alphaAPEsmear << "/";
    //   }

    //   // interpolator folder
    //   std::string subdir = "U_munu/";
    //   if (spatial_only) {
    //     subdir = "U_ij/";
    //   }
    //   oss << subdir;
    //   fsys::create_directories(fsys::absolute(oss.str()));

    //   oss << "C_glueball_";
    //   auto prevw = oss.width(8);
    //   auto prevf = oss.fill('0');
    //   oss << i;
    //   oss.width(prevw);
    //   oss.fill(prevf);

    //   const std::string path = oss.str();
    //   std::ofstream ofs(path, std::ios::out);

    //   std::vector<double> sinks[4]; // sink at 't'
    //   for (size_t i_PC = 0; i_PC < 4; i_PC++) {
    //     sinks[i_PC].resize(U.getLt());
    //   }

    //   for (size_t t = 0; t < U.getLt(); t++) {
    //     const accum Uij =
    //       operators::get_rest_tr_sum_U_munu<accum, Group>(U, t, spatial_only, false);
    //     const accum PUij =
    //       operators::get_rest_tr_sum_U_munu<accum, Group>(U, t, spatial_only, true);

    //     size_t i_PC = 0;
    //     for (int sP = 1; sP >= -1; sP -= 2) {
    //       for (int iC = 1; iC >= 0; iC--) {
    //         const bool C = bool(iC);
    //         std::complex<double> comb_snk = Uij + double(sP) * PUij;
    //         comb_snk /= 2.0;

    //         double snk;
    //         if (C) {
    //           snk = std::real(comb_snk);
    //         } else {
    //           snk = std::imag(comb_snk);
    //         }
    //         sinks[i_PC][t] = snk;

    //         ++i_PC;
    //       }
    //     }
    //   }

    //   ofs << "t C_{++}(t) phi_{++}(t) C_{+-}(t) phi_{+-}(t) C_{-+}(t) phi_{-+}(t) "
    //          "C_{--}(t) phi_{--}(t)"
    //       << std::endl; // header
    //   const size_t T_ext = U.getLt(); // lattice temporal time extent
    //   for (size_t t = 0; t < T_ext; t++) {
    //     ofs << t;
    //     for (size_t i_PC = 0; i_PC < 4; i_PC++) {
    //       double Ct = 0.0;
    //       for (size_t tau = 0; tau < T_ext; tau++) {
    //         Ct += sinks[i_PC][(t + tau) % T_ext] * sinks[i_PC][tau];
    //       }
    //       Ct /= double(T_ext); // average over all times
    //       ofs << " " << std::scientific << std::setprecision(16) << Ct << " "
    //           << sinks[i_PC][t];
    //     }
    //     ofs << std::endl;
    //   }

    //   ofs.close();

    //   return;
    // }

    // /**
    //  * @brief measure the glueball 0^{PC} correlators for the GEVP from a smeared gauge
    //  * configuration measure the glueball correlators ij for the 0^{++}, 0^{+-},
    //  0^{-+},
    //  * 0^{--} glueballs. The correlators are built as eq (1) of
    //  * https://arxiv.org/pdf/hep-lat/0603016.pdf
    //  * @tparam Group
    //  * @tparam sparams struct containing info on computation and output
    //  * @param U gauge configuration
    //  * @param i trajectory index
    //  * @param S specific parameters
    //  */
    // template <class Group, class sparams>
    // void meas_glueball_correlator_GEVP(const gaugeconfig<Group> &U,
    //                                    const size_t &i,
    //                                    const size_t &nAPEsmear,
    //                                    const sparams &S) {
    //   typedef typename accum_type<Group>::type accum;

    //   std::ostringstream oss_dir, oss_name;
    //   oss_dir << S.res_dir + "/";

    //   std::ostringstream oss_details;
    //   oss_details << "smearAPEn" << nAPEsmear << "alpha" << S.glueball.alphaAPEsmear;
    //   std::string meas_details = oss_details.str();

    //   if (S.glueball.doAPEsmear) {
    //     oss_dir << meas_details + "/";
    //   }
    //   fsys::create_directories(fsys::absolute(oss_dir.str()));

    //   if (S.glueball.doAPEsmear) {
    //     if (S.glueball.lengthy_file_name) {
    //       oss_name << "_" << meas_details;
    //     }
    //   }
    //   oss_name << "_";
    //   auto prevw = oss_name.width(8);
    //   auto prevf = oss_name.fill('0');
    //   oss_name << i;
    //   oss_name.width(prevw);
    //   oss_name.fill(prevf);

    //   const size_t T_ext = U.getLt(); // lattice temporal time extent
    //   const size_t rmin = S.glueball.rmin_GEVP, rmax = S.glueball.rmax_GEVP;
    //   const size_t N_ops = rmax; // number of GEVP interpolatoing operators

    //   // phi_i(t)^{PC}
    //   std::vector<std::vector<std::array<std::array<double, 2>, 2>>> phi(rmax - rmin +
    //   1);

    //   for (size_t ir = 0; ir <= rmax - rmin; ir++) {
    //     const size_t r = ir + rmin;
    //     phi[ir].resize(T_ext);
    //     for (size_t t = 0; t < T_ext; t++) {
    //       const std::complex<double> Pp =
    //         glueballs::rest_trace_wloop_munu<double, Group>(t, U, r, r, 1, 2, false);
    //       const std::complex<double> Pm =
    //         glueballs::rest_trace_wloop_munu<double, Group>(t, U, r, r, 1, 2, true);
    //       phi[ir][t][0][0] = (Pp + Pm).real() / 2.0; // PC=++
    //       phi[ir][t][0][1] = (Pp + Pm).imag() / 2.0; // PC=+-
    //       phi[ir][t][1][0] = (Pp - Pm).real() / 2.0; // PC=-+
    //       phi[ir][t][1][1] = (Pp - Pm).imag() / 2.0; // PC=--
    //     }
    //   }

    //   for (size_t i1 = 0; i1 <= rmax - rmin; i1++) {
    //     const size_t r1 = i1 + rmin;
    //     for (size_t i2 = 0; i2 <= i1; i2++) { // C_{ij} == C_{ji}
    //       const size_t r2 = i2 + rmin;
    //       const std::string dir_ij = oss_dir.str() + std::to_string(r1) + "_" +
    //                                  std::to_string(r2) + "/"; // directory path
    //       fsys::create_directories(fsys::absolute(dir_ij)); // creating directory

    //       const std::string path =
    //         dir_ij + "C_glueball" + oss_name.str(); // full path of output file

    //       std::ostringstream oss_ij;
    //       oss_ij << "t "
    //              << "C_{++}(t) phi_i_{++}(t) phi_j_{++}(t) "
    //              << "C_{+-}(t) phi_i_{+-}(t) phi_j_{+-}(t) "
    //              << "C_{-+}(t) phi_i_{-+}(t) phi_j_{-+}(t) "
    //              << "C_{--}(t) phi_i_{--}(t) phi_j_{--}(t)" << std::endl; // header

    //       for (size_t t = 0; t < T_ext; t++) {
    //         oss_ij << t;
    //         for (size_t P = 0; P <= 1; P++) {
    //           for (size_t C = 0; C <= 1; C++) {
    //             double Ct = 0.0;
    //             for (size_t tau = 0; tau < T_ext; tau++) {
    //               Ct += phi[i1][(t + tau) % T_ext][P][C] * phi[i2][tau][P][C];
    //             }
    //             Ct /= double(T_ext); // average over all times
    //             oss_ij << " " << std::scientific << std::setprecision(16) << Ct << " "
    //                    << phi[i1][t][P][C] << " " << phi[i2][t][P][C];
    //           }
    //         }
    //         oss_ij << std::endl;
    //       }

    //       std::ofstream ofs(path, std::ios::out);
    //       ofs << oss_ij.str();
    //       ofs.close();
    //     }
    //   }

    //   return;
    // }

    /**
     * @brief measure the glueball 0^{PC} interpolators (at rest) from a smeared gauge
     * configuration
     * @tparam Group
     * @tparam sparams struct containing info on computation and output
     * @param type type of interpolating field
     * @param U gauge configuration
     * @param i trajectory index
     * @param S specific parameters
     */
    template <class Group, class sparams>
    void meas_glueball_interpolators(const std::string &type,
                                     const gaugeconfig<Group> &U,
                                     const size_t &i,
                                     const size_t &nAPEsmear,
                                     const sparams &S) {
      typedef typename accum_type<Group>::type accum;

      std::ostringstream oss_dir, oss_name;
      oss_dir << S.res_dir + "/glueball/interpolator/" + type + "/";

      std::ostringstream oss_details;
      oss_details << "smearAPEn" << nAPEsmear << "alpha" << S.glueball.alphaAPEsmear;
      std::string meas_details = oss_details.str();

      if (S.glueball.doAPEsmear) {
        oss_dir << meas_details + "/";
      }
      fsys::create_directories(fsys::absolute(oss_dir.str()));

      if (S.glueball.doAPEsmear) {
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

      const size_t T_ext = U.getLt(); // lattice temporal time extent
      const size_t rmin = S.glueball.rmin, rmax = S.glueball.rmax;
      const size_t N_ops = rmax; // number of GEVP interpolatoing operators

      // phi_i(t)^{PC}
      std::vector<std::vector<std::array<std::array<double, 2>, 2>>> phi(rmax - rmin + 1);

      for (size_t i1 = 0; i1 <= rmax - rmin; i1++) {
        const size_t r1 = i1 + rmin;
        for (size_t i2 = 0; i2 <= i1; i2++) { // phi_{ij} == phi_{ji}
          const size_t r2 = i2 + rmin;

          const std::string dir_ij = oss_dir.str() + std::to_string(r1) + "_" +
                                     std::to_string(r2) + "/"; // directory path
          fsys::create_directories(fsys::absolute(dir_ij)); // creating directory

          const std::string path =
            dir_ij + "phi" + oss_name.str(); // full path of output file

          std::ostringstream oss_ij;
          oss_ij << "t pp pm mp mm" << std::endl; // header

          for (size_t t = 0; t < T_ext; t++) {
            oss_ij << t;
            const std::complex<double> Pp = glueballs::get_rest_trace_loop<double, Group>(
              type, t, U, r1, r2, false, S.glueball.spatial_loops);
            const std::complex<double> Pm = glueballs::get_rest_trace_loop<double, Group>(
              type, t, U, r1, r2, true, S.glueball.spatial_loops);
            double phi[2][2];
            phi[0][0] = (Pp + Pm).real() / 2.0; // PC=++
            phi[0][1] = (Pp + Pm).imag() / 2.0; // PC=+-
            phi[1][0] = (Pp - Pm).real() / 2.0; // PC=-+
            phi[1][1] = (Pp - Pm).imag() / 2.0; // PC=--
            for (size_t P = 0; P <= 1; P++) {
              for (size_t C = 0; C <= 1; C++) {
                oss_ij << " " << std::scientific << std::setprecision(16) << phi[P][C];
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


    // /**
    //  * @brief measure the glueball 0^{PC} correlators for the GEVP from a smeared gauge configuration measure the glueball correlators ij for the 0^{++}, 0^{+-}, 0^{-+}, 0^{--} glueballs. The correlators are built as eq (1) of
    //  * https://arxiv.org/pdf/hep-lat/0603016.pdf
    //  * @tparam Group
    //  * @tparam sparams struct containing info on computation and output
    //  * @param U gauge configuration
    //  * @param i trajectory index
    //  * @param S specific parameters
    //  */
    // template <class Group, class sparams>
    // void meas_glueball_correlator(const std::string& type,const gaugeconfig<Group> &U,
    //                                    const size_t &i,
    //                                    const size_t &nAPEsmear,
    //                                    const sparams &S) {
    //   std::ostringstream oss_dir, oss_name;
    //   oss_dir << S.res_dir + "/glueball/correlator/";

    //   std::ostringstream oss_details;
    //   oss_details << "smearAPEn" << nAPEsmear << "alpha" << S.glueball.alphaAPEsmear;
    //   std::string meas_details = oss_details.str();

    //   if (S.glueball.doAPEsmear) {
    //     oss_dir << meas_details + "/";
    //   }
    //   fsys::create_directories(fsys::absolute(oss_dir.str()));

    //   if (S.glueball.doAPEsmear) {
    //     if (S.glueball.lengthy_file_name) {
    //       oss_name << "_" << meas_details;
    //     }
    //   }
    //   oss_name << "_";
    //   auto prevw = oss_name.width(8);
    //   auto prevf = oss_name.fill('0');
    //   oss_name << i;
    //   oss_name.width(prevw);
    //   oss_name.fill(prevf);

    //   const size_t T_ext = U.getLt(); // lattice temporal time extent
    //   const size_t rmin = S.glueball.rmin_GEVP, rmax = S.glueball.rmax_GEVP;
    //   const size_t N_ops = rmax; // number of GEVP interpolatoing operators

    //   // phi_i(t)^{PC}
    //   std::vector<std::vector<std::array<std::array<double, 2>, 2>>> phi(rmax - rmin +
    //   1);

    //   for (size_t ir = 0; ir <= rmax - rmin; ir++) {
    //     const size_t r = ir + rmin;
    //     phi[ir].resize(T_ext);
    //     for (size_t t = 0; t < T_ext; t++) {
    //       const std::complex<double> Pp =
    //         glueballs::rest_trace_wloop_munu<double, Group>(t, U, r, r, 1, 2, false);
    //       const std::complex<double> Pm =
    //         glueballs::rest_trace_wloop_munu<double, Group>(t, U, r, r, 1, 2, true);
    //       phi[ir][t][0][0] = (Pp + Pm).real() / 2.0; // PC=++
    //       phi[ir][t][0][1] = (Pp + Pm).imag() / 2.0; // PC=+-
    //       phi[ir][t][1][0] = (Pp - Pm).real() / 2.0; // PC=-+
    //       phi[ir][t][1][1] = (Pp - Pm).imag() / 2.0; // PC=--
    //     }
    //   }

    //   for (size_t i1 = 0; i1 <= rmax - rmin; i1++) {
    //     const size_t r1 = i1 + rmin;
    //     for (size_t i2 = 0; i2 <= i1; i2++) { // C_{ij} == C_{ji}
    //       const size_t r2 = i2 + rmin;
    //       const std::string dir_ij = oss_dir.str() + std::to_string(r1) + "_" +
    //                                  std::to_string(r2) + "/"; // directory path
    //       fsys::create_directories(fsys::absolute(dir_ij)); // creating directory

    //       const std::string path =
    //         dir_ij + "C_glueball" + oss_name.str(); // full path of output file

    //       std::ostringstream oss_ij;
    //       oss_ij << "t "
    //              << "C_{++}(t) phi_i_{++}(t) phi_j_{++}(t) "
    //              << "C_{+-}(t) phi_i_{+-}(t) phi_j_{+-}(t) "
    //              << "C_{-+}(t) phi_i_{-+}(t) phi_j_{-+}(t) "
    //              << "C_{--}(t) phi_i_{--}(t) phi_j_{--}(t)" << std::endl; // header

    //       for (size_t t = 0; t < T_ext; t++) {
    //         oss_ij << t;
    //         for (size_t P = 0; P <= 1; P++) {
    //           for (size_t C = 0; C <= 1; C++) {
    //             double Ct = 0.0;
    //             for (size_t tau = 0; tau < T_ext; tau++) {
    //               Ct += phi[i1][(t + tau) % T_ext][P][C] * phi[i2][tau][P][C];
    //             }
    //             Ct /= double(T_ext); // average over all times
    //             oss_ij << " " << std::scientific << std::setprecision(16) << Ct << " "
    //                    << phi[i1][t][P][C] << " " << phi[i2][t][P][C];
    //           }
    //         }
    //         oss_ij << std::endl;
    //       }

    //       std::ofstream ofs(path, std::ios::out);
    //       ofs << oss_ij.str();
    //       ofs.close();
    //     }
    //   }

    //   return;
    // }


  } // namespace from_smeared_field

  // /**
  //  * @brief measure the glueball 0^{PC} correlators for all numbers of smearing steps
  //  *
  //  * @tparam Group
  //  * @tparam sparams struct containing info on computation and output
  //  * @param U gauge configuration
  //  * @param i trajectory index
  //  * @param S specific parameters
  //  */
  // template <class Group, class sparams>
  // void meas_glueball_correlator_U_munu(const gaugeconfig<Group> &U0,
  //                                      const size_t &i,
  //                                      const sparams &S,
  //                                      const bool &spatial_only) {
  //   const std::vector<size_t> &v_ns =
  //     S.glueball.vec_nAPEsmear; // vector of number of smearing steps
  //   const size_t ns = v_ns.size();

  //   gaugeconfig<Group> U = U0; // copy of the initial gauge configuration
  //   for (size_t is = 0; is < ns; is++) {
  //     const size_t nsteps = (is == 0) ? v_ns[is] : v_ns[is] - v_ns[is - 1];
  //     for (size_t i = 0; i < nsteps; i++) {
  //       // other nsteps so that we reach the value of v_ns[is]
  //       spatial_APEsmearing<double, Group>(U, S.glueball.alphaAPEsmear);
  //     }
  //     from_smeared_field::meas_glueball_correlator_U_munu(U, i, v_ns[is], S,
  //                                                         spatial_only);
  //   }
  // }

  // /**
  //  * @brief measure of the glueball correlators for all numbers of smearing steps
  //  *
  //  * @tparam Group
  //  * @tparam sparams
  //  * @param U0 initial gauge configuration
  //  * @param i trajectory index
  //  * @param S specific parameters
  //  */
  // template <class Group, class sparams>
  // void meas_glueball_correlator_GEVP(const gaugeconfig<Group> &U0,
  //                                    const size_t &i,
  //                                    const sparams &S) {
  //   const std::vector<size_t> &v_ns =
  //     S.glueball.vec_nAPEsmear; // vector of number of smearing steps
  //   const size_t ns = v_ns.size();

  //   gaugeconfig<Group> U = U0; // copy of the initial gauge configuration
  //   for (size_t is = 0; is < ns; is++) {
  //     const size_t nsteps = (is == 0) ? v_ns[is] : v_ns[is] - v_ns[is - 1];
  //     for (size_t i = 0; i < nsteps; i++) {
  //       // other nsteps so that we reach the value of v_ns[is]
  //       spatial_APEsmearing<double, Group>(U, S.glueball.alphaAPEsmear);
  //     }
  //     from_smeared_field::meas_glueball_correlator_GEVP(U, i, v_ns[is], S);
  //   }
  // }

  /**
   * @brief measure of the glueball interpolators for all numbers of smearing
   * steps
   *
   * @tparam Group
   * @tparam sparams
   * @param type type of interpolator
   * @param U0 initial gauge configuration
   * @param i trajectory index
   * @param S specific parameters
   */
  template <class Group, class sparams>
  void meas_glueball_interpolators(const std::string &type,
                                   const gaugeconfig<Group> &U0,
                                   const size_t &i,
                                   const sparams &S) {
    const std::vector<size_t> &v_ns =
      S.glueball.vec_nAPEsmear; // vector of number of smearing steps
    const size_t ns = v_ns.size();

    gaugeconfig<Group> U = U0; // copy of the initial gauge configuration
    for (size_t is = 0; is < ns; is++) {
      const size_t nsteps = (is == 0) ? v_ns[is] : v_ns[is] - v_ns[is - 1];
      for (size_t i = 0; i < nsteps; i++) {
        // other nsteps so that we reach the value of v_ns[is]
        spatial_APEsmearing<double, Group>(U, S.glueball.alphaAPEsmear);
      }
      from_smeared_field::meas_glueball_interpolators(type, U, i, v_ns[is], S);
    }
  }

  /**
   * @brief measure of the glueball correlator
   *
   * @tparam Group
   * @tparam sparams
   * @param type type of interpolator used
   * @param U0 initial gauge configuration
   * @param i trajectory index
   * @param S specific parameters
   */
  template <class Group, class sparams>
  void meas_glueball_correlator(const std::string&  type,const gaugeconfig<Group> &U0,
                                     const size_t &i,
                                     const sparams &S) {
    const std::vector<size_t> &v_ns =
      S.glueball.vec_nAPEsmear; // vector of number of smearing steps
    const size_t ns = v_ns.size();

    gaugeconfig<Group> U = U0; // copy of the initial gauge configuration
    for (size_t is = 0; is < ns; is++) {
      const size_t nsteps = (is == 0) ? v_ns[is] : v_ns[is] - v_ns[is - 1];
      for (size_t i = 0; i < nsteps; i++) {
        // other nsteps so that we reach the value of v_ns[is]
        spatial_APEsmearing<double, Group>(U, S.glueball.alphaAPEsmear);
      }
      from_smeared_field::meas_glueball_correlator(type, U, i, v_ns[is], S);
    }
  }


  /**
   * measures the planar wilson loops in temporal and spatial direction
   * -(0,1) loops for temporal, (1,2) and (1,3) loops for spatial
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
                             const double &sizeWloops,
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
   * measures the nonplanar wilson loops in temporal and spatial direction
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
                                const double &sizeWloops,
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
