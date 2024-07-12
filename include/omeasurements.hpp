// omeasurements.hpp
/**
 * @brief routines for online/offline measurements of various observables.
 * This file defines functions which handle the measure and printing of observables over a
 * given gauge configuration.
 */

#pragma once
#include <complex>

#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "glueballs.hpp"
#include "gradient_flow.hh"
#include "io.hh"
#include "links.hpp"
#include "obc_gaugemonomial.hh"
#include "operators.hpp"
#include "parameters.hh"
#include "propagator.hpp"
#include "smearape.hh"
#include "wilsonloop.hh"
#include "polyakov_loop.hh"

#include <boost/filesystem.hpp>

#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
using Complex = std::complex<double>;

namespace omeasurements {

  namespace fsys = boost::filesystem;

  /**
   * @brief Get the retr realtrace density object
   * 
   * @tparam Gauge group
   * @param U gauge configuration
   * @return double 
   */
  template <class Group>
  double get_retr_realtrace_density(const gaugeconfig<Group> &U, const std::string &bc){
  //double denuminator = U.getVolume()*U.getndims();
  //std::cout << "gaugemass in retrace" << gaugemass << "\n ";
  obc::weights w(bc, U.getLx(), U.getLy(), U.getLz(), U.getLt(), U.getndims());
  double realtrace = obc::retr_sum_realtrace(U, w);
  double num_lattice_sites = 1;
  const size_t ndims = U.getndims();
  if (ndims > 0) {
    num_lattice_sites *= U.getLt();
    if (ndims > 1) {
      num_lattice_sites *= U.getLx();
      if (ndims > 2){
        num_lattice_sites *= U.getLy();
        if (ndims > 3){
          num_lattice_sites *= U.getLz();
        }
      }
    }
  }
  double normalized_retrace =  (1./(2. * num_lattice_sites*ndims)) * realtrace; 
  return normalized_retrace;
  //realtrace /= denuminator;
  //return realtrace;
  }

  /**
   * @brief Get the retr plaquette density object
   *
   * @param bc
   * @return double
   */
  template <class Group>
  double get_retr_plaquette_density(const gaugeconfig<Group> &U,
                                    const std::string &bc,
                                    const bool &spatial = false) {
    const double ndims_fact =
      spacetime_lattice::num_pLloops_half(U.getndims() - int(spatial));
    double P = 0.0; // plaquette density value
    double den = ndims_fact * double(U.getNc());

    if (bc == "periodic") {
      P = flat_spacetime::retr_sum_Wplaquettes(U, /*xi=*/1.0, /*anisotropic=*/false,
                                               /*spatial=*/spatial);
      den *= U.getVolume();
    } else if (bc == "spatial_open") {
      const size_t ndims = U.getndims();
      obc::weights w(bc, U.getLx(), U.getLy(), U.getLz(), U.getLt(), ndims);

      P = obc::retr_sum_Wplaquettes(U, w, /*xi=*/1.0, /*anisotropic=*/false,
                                    /*spatial=*/spatial);
      if (ndims > 0) {
        den *= U.getLt();
        if (ndims > 1) {
          den *= (U.getLx() - 1);
          if (ndims > 2) {
            den *= (U.getLy() - 1);
            if (ndims > 3) {
              den *= (U.getLz() - 1);
            }
          }
        }
      }
    }

    P /= den; // normalizing
    return P;
  }

/**
 * @brief calculates the realtrace density and saves it in a file 
 * 
 * @tparam Group 
 * @tparam sparams struct containing info on the computation and physics
 * @param U gauge configuration
 * @param i configuration index 
 * @param pparams physical parameters
 * @param S realtrace configuration paramters
 */
template <class Group, class sparams>
void meas_realtrace(const gaugeconfig <Group> U,
                    const size_t &i,
                    const global_parameters::physics &pparams,
                    const sparams &S){
  std::ostringstream oss;
  oss << S.res_dir +"/" + S.retrace.subdir <<"/";
  fsys::create_directories(fsys::absolute(oss.str()));

  oss << "retrace_";
  auto prevw = oss.width(8);
  auto prevf = oss.fill('0');
  oss << i;
  oss.width(prevw);
  oss.fill(prevf);

  const std::string path = oss.str();
  std::ofstream ofs(path, std::ios::out);

  ofs << "i retrace \n";

  const double retrace = get_retr_realtrace_density(U, S.retrace.bc);
  ofs << i << std::scientific << std::setprecision(16) << " " << retrace << "\n";
  ofs.close();

  return;
}

  template <class Group, class sparams>
  void meas_plaquette(const gaugeconfig<Group> &U,
                      const size_t &i,
                      const global_parameters::physics &pparams,
                      const sparams &S) {
    std::ostringstream oss;
    oss << S.res_dir + "/" + S.plaquette.subdir << "/";
    fsys::create_directories(fsys::absolute(oss.str())); // creating directory

    oss << "plaquette_";
    auto prevw = oss.width(8);
    auto prevf = oss.fill('0');
    oss << i;
    oss.width(prevw);
    oss.fill(prevf);

    const std::string path = oss.str();
    std::ofstream ofs(path, std::ios::out);

    ofs << "i P\n";
    // plaquette density value
    const double P = get_retr_plaquette_density(U, S.plaquette.bc, S.plaquette.spatial);
    ofs << i << std::scientific << std::setprecision(16) << " " << P << "\n";
    ofs.close();

    return;
  }

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
  template <class Group, class sparams>
  void meas_spatial_polyakov(const gaugeconfig<Group> &U,
                                          const size_t &i,
                                          const sparams &S){
/*  std::ostringstream os;
  os << res_dir + "/polyakovloop.";

  auto prevw = os.width(6);
  auto prevf = os.fill('0');
  os << i;
  os.width(prevw);
  os.fill(prevf);

*/ 
std::ostringstream oss;
oss << S.res_dir +"/" + S.spatial_polyakov.subdir <<"/";
fsys::create_directories(fsys::absolute(oss.str()));

oss << "polyakov_";
auto prevw = oss.width(8);
auto prevf = oss.fill('0');
oss << i;
oss.width(prevw);
oss.fill(prevf);

const std::string path = oss.str();
std::ofstream ofs(path, std::ios::out);

ofs << "i polyreal polyimag \n";

const Complex spatial_polyakov = polyakov_loop::polyakov_loop_spatial_average(U);
ofs << i << std::scientific << std::setprecision(16) << " " << std::real(spatial_polyakov) << " " << std::imag(spatial_polyakov) << "\n";
ofs.close();
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
    const std::string res_dir = S.res_dir + "/" + S.gradient_flow.subdir + "/";
    fsys::create_directories(fsys::absolute(res_dir));

    const double eps = S.gradient_flow.epsilon;
    const double tmax = S.gradient_flow.tmax;
    double tstart = S.gradient_flow.tstart;
    const bool save_conf = S.gradient_flow.save_conf;

    std::ostringstream os;
    os << res_dir + "/gradient_flow.";
    auto prevw = os.width(6);
    auto prevf = os.fill('0');
    os << i;
    os.width(prevw);
    os.fill(prevf);

    gaugeconfig<Group> V = U;
    if (tstart > eps) {
      V.load(os.str() + "_t" + std::to_string(tstart) + ".conf");
    }
    flat_spacetime::gradient_flow(V, os.str(), tmax, eps, pparams.xi, tstart, save_conf);

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

    /**
     * @brief measure the glueball 0^{PC} interpolators (at rest) from a smeared gauge
     * configuration
     * @tparam Group
     * @tparam sparams struct containing info on computation and output
     * @param type type of interpolating field
     * @param U gauge configuration
     * @param i trajectory index
     * @param dump true when the interpolator values are printed out
     * @param S specific parameters
     */
    template <class Group, class sparams>
    std::vector<xt::xarray<double>>
    meas_glueball_interpolators(const std::string &type,
                                const gaugeconfig<Group> &U,
                                const size_t &i,
                                const size_t &nAPEsmear,
                                const bool &save_interpolator,
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
      //      fsys::create_directories(fsys::absolute(oss_dir.str()));

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
      const size_t N_ops = rmax; // number of interpolatoing operators

      // phi_i(t)^{PC}
      const size_t nr = (rmax - rmin + 1);
      std::vector<xt::xarray<double>> phi(nr); // vector of all interpolators
      size_t ii = 0; // interpolator index

      for (size_t i1 = 0; i1 < nr; i1++) {
        const size_t r1 = i1 + rmin;

        phi[ii].resize({T_ext, 5}); // t, ++, +-, -+, --

        for (size_t t = 0; t < T_ext; t++) {
          phi[ii](t, 0) = t;
          const std::complex<double> Pp = glueballs::get_rest_trace_loop<double, Group>(
            type, t, U, r1, false, S.glueball.spatial_loops);
          const std::complex<double> Pm = glueballs::get_rest_trace_loop<double, Group>(
            type, t, U, r1, true, S.glueball.spatial_loops);
          phi[ii](t, 1) = (Pp + Pm).real() / 2.0; // PC=++
          phi[ii](t, 3) = (Pp + Pm).imag() / 2.0; // PC=+-
          phi[ii](t, 2) = (Pp - Pm).real() / 2.0; // PC=-+
          phi[ii](t, 4) = (Pp - Pm).imag() / 2.0; // PC=--
        }

        if (save_interpolator) {
          // directory path
          const std::string dir_ij = oss_dir.str() + std::to_string(r1) + "/";
          fsys::create_directories(fsys::absolute(dir_ij));
          // full path of output file
          const std::string path = dir_ij + "phi" + oss_name.str();
          std::ofstream ofs(path, std::ios::out);
          ofs << std::scientific << std::setprecision(16);
          io::xtensor_to_stream(ofs, phi[ii], " ", "t pp pm mp mm");
          ofs.close();
        }
        ++ii;
      }

      return phi;
    }

    /**
     * @brief measure the glueball 0^{PC} correlators
     * measure the glueball correlators ij for the 0^{++}, 0^{--} glueballs.
     * @tparam Group
     * @tparam sparams struct containing info on computation and output
     * @param U gauge configuration
     * @param i trajectory index
     * @param S specific parameters
     */
    template <class Group, class sparams>
    void meas_glueball_correlator(const std::string &type,
                                  const gaugeconfig<Group> &U,
                                  const size_t &i,
                                  const size_t &nAPEsmear,
                                  const sparams &S) {
      std::ostringstream oss_dir, oss_name;
      oss_dir << S.res_dir + "/glueball/correlator/" + type + "/";

      std::ostringstream oss_details;
      oss_details << "smearAPEn" << nAPEsmear << "alpha" << S.glueball.alphaAPEsmear;
      std::string meas_details = oss_details.str();

      if (S.glueball.doAPEsmear) {
        oss_dir << meas_details + "/";
      }
      //      fsys::create_directories(fsys::absolute(oss_dir.str()));

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
      const size_t N_ops = rmax; // number of interpolatoing operators

      // phi_i(t)^{PC}
      std::vector<xt::xarray<double>> phi = meas_glueball_interpolators(
        type, U, i, nAPEsmear, S.glueball.save_interpolator, S);
      if (!S.glueball.correlator) {
        return;
      }

      const size_t n_phi = phi.size();

      for (size_t i = 0; i < n_phi; i++) {
        for (size_t j = 0; j <= i; j++) {
          // directory path
          const std::string dir_ij = oss_dir.str() + std::to_string(i + S.glueball.rmin) +
                                     "_" + std::to_string(j + S.glueball.rmin) + "/";
          fsys::create_directories(fsys::absolute(dir_ij)); // creating directory

          // full path of output file
          const std::string path = dir_ij + "C_glueball" + oss_name.str();

          const std::string header = "t pp pm mp mm";
          xt::xarray<double> corr;
          corr.resize({T_ext, 5}); // t, ++, +-, -+, --

          for (size_t t = 0; t < T_ext; t++) {
            corr(t, 0) = t;
            size_t iPC = 1; // 1st column is time 't'
            for (size_t P = 0; P <= 1; P++) {
              for (size_t C = 0; C <= 1; C++) {
                corr(t, iPC) = 0.0;
                for (size_t tau = 0; tau < T_ext; tau++) {
                  const size_t t1 = (t + tau) % T_ext;
                  corr(t, iPC) += phi[i](t1, iPC) * phi[j](tau, iPC);
                }
                corr(t, iPC) /= double(T_ext); // average over all times
                ++iPC;
              }
            }
          }

          std::ofstream ofs(path, std::ios::out);
          ofs << std::scientific << std::setprecision(16);
          io::xtensor_to_stream(ofs, corr, " ", header);
          ofs.close();
        }
      }

      return;
    }

  } // namespace from_smeared_field

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
      auto foo = from_smeared_field::meas_glueball_interpolators(
        type, U, i, v_ns[is], S.glueball.save_interpolator, S);
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
  void meas_glueball_correlator(const std::string &type,
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
 * @brief saves the correlators of the polyakov loops to file
 * 
 * @param U pointer to config
 * @param S pointer to measure params
 * @param i index of configuration 
 */
 template <class Group, class sparams>
  void meas_correlator_polyakov_pot(const gaugeconfig<Group> &U,
                                    const sparams &S,
                                    const size_t &i
                                    ){
                                      const size_t rmax = S.polypot.rmax;
                                      Complex correlator;
                                      std::ostringstream oss;
                                      oss << S.res_dir + "/" + S.polypot.subdir +"/";
                                      fsys::create_directories(fsys::absolute(oss.str()));

                                      oss << "polypot_";
                                      auto prevw = oss.width(8);
                                      auto prevf = oss.fill('0');
                                      oss << i;
                                      oss.width(prevw);
                                      oss.fill(prevf);

                                      const std::string path = oss.str();
                                      std::ofstream ofs(path, std::ios::out);

                                      ofs << "i ";
                                      for(size_t r= 1; r <= rmax; r++){
                                        ofs << std::to_string(r) << " ";
                                      }
                                      ofs << "\n";
                                      ofs << i << " ";
                                      //resultfile.open(filename_polyakov_corr, std::ios::app);
                                      //resultfile << i;
                                      for(size_t r = 1; r <= rmax; r++){
                                        correlator = polyakov_loop::polyakov_loop_correlator(U, r);
                                        ofs << std::scientific << std::setprecision(16) << std::real(correlator) << " ";
                                        //resultfile << std::setw(14) << std::scientific << correlator << " ";
                                      }
                                      ofs << "\n";
                                      ofs.close();
                                      return;
                                      //resultfile << std::endl;
                                      //resultfile.close() 
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
