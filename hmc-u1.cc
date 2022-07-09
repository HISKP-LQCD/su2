/**
 * @file hmc-u1.cc
 * @author Carsten Urbach (urbach@hiskp.uni-bonn.de)
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief Hybrid Monte Carlo for a U(1) theory
 * @version 0.1
 * @date 2022-05-11
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "flat-energy_density.hh"
#include "flat-gauge_energy.hpp"
#include "gaugeconfig.hh"
#include "integrator.hh"
#include "io.hh"
#include "md_update.hh"
#include "monomial.hh"
#include "omeasurements.hpp"
#include "parse_input_file.hh"
#include "random_gauge_trafo.hh"
#include "su2.hh"
#include "u1.hh"
#include "version.hh"

#include "rotating-gaugemonomial.hpp"

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

#include "detDDdag_monomial.hh"

namespace po = boost::program_options;
namespace gp = global_parameters;
namespace fsys = boost::filesystem;
/**
 * @brief initialize gaige configuration for the hmc
 */
void initialize_U_hmc(gaugeconfig<_u1> &U,
                      bool &g_heat,
                      size_t &g_icounter,
                      double &normalisation,
                      const gp::physics &pparams,
                      const gp::hmc_u1 &hparams);

int main(int ac, char *av[]) {
  std::cout << "## HMC Algorithm for U(1) gauge theory\n";
  std::cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017, 2021)\n";
  std::cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << "\n";

  gp::physics pparams; // physics parameters
  gp::hmc_u1 hparams; // hmc parameters

  std::string input_file; // yaml input file path
  input_file_parsing::parse_command_line(ac, av, input_file);

  namespace in_hmc = input_file_parsing::u1::hmc;
  in_hmc::parse_input_file(input_file, pparams, hparams);

  // configurations folder
  fsys::create_directories(fsys::absolute(hparams.conf_dir));

  // online measurements folder
  const fsys::path abs_omeas_res = fsys::absolute(hparams.omeas.res_dir);
  fsys::create_directories(abs_omeas_res);

  bool g_heat; // hot or cold starting configuration
  size_t g_icounter; // 1st configuration(trajectory) to load from
  gaugeconfig<_u1> U(pparams.Lx, pparams.Ly, pparams.Lz, pparams.Lt, pparams.ndims,
                     pparams.beta);
  double normalisation;
  initialize_U_hmc(U, g_heat, g_icounter, normalisation, pparams, hparams);

  // Molecular Dynamics parameters
  md_params mdparams(hparams.n_steps, hparams.tau);

  // generate list of monomials
  std::list<monomial<double, _u1> *> monomial_list;
  flat_spacetime::gaugemonomial<double, _u1> gm(0, pparams.xi);
  rotating_spacetime::gauge_monomial<double, _u1> gm_rot(0, pparams.Omega);

  kineticmonomial<double, _u1> km(0);
  km.setmdpassive();
  monomial_list.push_back(&km);

  staggered::detDDdag_monomial<double, _u1> detDDdag(
    0, pparams.m0, hparams.solver, hparams.tolerance_cg, hparams.seed_pf,
    hparams.solver_verbosity);

  if (pparams.include_gauge) {
    if (pparams.rotating_frame) {
      monomial_list.push_back(&gm_rot);
    } else {
      monomial_list.push_back(&gm);
    }
  }

  if (pparams.include_staggered_fermions) { // including S_F (fermionic) in the action
    monomial_list.push_back(&detDDdag);
  }

  // setting up the integrator
  integrator<double, _u1> *md_integ =
    set_integrator<double, _u1>(hparams.integrator, hparams.exponent);

  std::ofstream os;
  if (g_icounter == 0) {
    os.open(hparams.conf_dir + "/output.hmc.data", std::ios::out);
  } else {
    os.open(hparams.conf_dir + "/output.hmc.data", std::ios::app);
  }
  std::cout << "## Normalization factor: A = 2/(d*(d-1)*N_lat*N_c) = " << std::scientific
            << std::setw(18) << std::setprecision(15) << normalisation << "\n";
  std::cout << "## Acceptance rate parcentage: rho = rate/(i+1)\n";

  // header: column names in the output
  std::string head_str = io::hmc::get_header(" ");
  std::cout << head_str;
  os << head_str;

  double rate = 0.;

  const std::string conf_path_basename = io::get_conf_path_basename(pparams, hparams);

  for (size_t i = g_icounter; i < hparams.n_meas + g_icounter; i++) {
    if (hparams.do_hmc) {
      mdparams.disablerevtest();
      if (i > 0 && hparams.N_rev != 0 && (i) % hparams.N_rev == 0) {
        mdparams.enablerevtest();
      }
      // PRNG engine
      std::mt19937 engine(hparams.seed + i);
      // perform the MD update
      md_update(U, engine, mdparams, monomial_list, *md_integ);

      const double energy = flat_spacetime::gauge_energy(U);
      double E = 0., Q = 0.;
      flat_spacetime::energy_density(U, E, Q);
      rate += mdparams.getaccept();

      std::cout << i << " " << mdparams.getaccept() << " " << std::scientific
                << std::setw(18) << std::setprecision(15) << energy * normalisation << " "
                << std::setw(15) << mdparams.getdeltaH() << " " << std::setw(15)
                << rate / static_cast<double>(i + 1) << " ";

      if (mdparams.getrevtest()) {
        std::cout << mdparams.getdeltadeltaH();
      } else {
        std::cout << "NA";
      }
      std::cout << " " << Q << std::endl;

      os << i << " " << mdparams.getaccept() << " " << std::scientific << std::setw(18)
         << std::setprecision(15) << energy * normalisation << " " << std::setw(15)
         << mdparams.getdeltaH() << " " << std::setw(15)
         << rate / static_cast<double>(i + 1) << " ";
      if (mdparams.getrevtest()) {
        os << mdparams.getdeltadeltaH();
      } else {
        os << "NA";
      }
      os << " " << Q << std::endl;
    }

    if (i > 0 && (i % hparams.N_save) == 0) { // saving U after each N_save trajectories
      std::string path_i = conf_path_basename + "." + std::to_string(i);

      if (hparams.do_hmc) {
        U.save(path_i);
      } else {
        int lerr = U.load(path_i);
        if (lerr == 1) {
          continue;
        }
      }

      // online measurements
      bool do_omeas =
        hparams.make_omeas && i > hparams.omeas.icounter && i % hparams.omeas.nstep == 0;
      if (hparams.do_hmc) {
        // check also if trajectory was accepted
        do_omeas = do_omeas && mdparams.getaccept();
      }

      if (do_omeas) {
        if (i == g_icounter && hparams.do_hmc) {
          continue; // online measurements already done
        }

        if (hparams.omeas.Wloop) {
          if (hparams.omeas.verbosity > 0) {
            std::cout << "## online measuring: Wilson loop\n";
          }
          omeasurements::meas_wilson_loop<_u1>(U, i, hparams.conf_dir);
        }
        if (hparams.omeas.gradient_flow) {
          if (hparams.omeas.verbosity > 0) {
            std::cout << "## online measuring: Gradient flow\n";
          }
          omeasurements::meas_gradient_flow<_u1>(U, i, hparams.omeas);
        }

        if (hparams.omeas.pion_staggered) {
          if (hparams.omeas.verbosity > 0) {
            std::cout << "## online measuring: Pion correlator\n";
          }
          omeasurements::meas_pion_correlator<_u1>(U, i, pparams.m0, hparams.omeas);
        }

        if (hparams.omeas.measure_glueball_params.do_measure) {
          if (hparams.omeas.verbosity > 0) {
            std::cout << "## online measuring: Glueballs 0^{PC} correlators\n";
          }
          omeasurements::meas_glueball_correlator<_u1>(U, i, hparams.omeas);
        }
      }

      if (hparams.do_hmc) { // storing last conf index (only after online measurements has
                            // been done)
        io::hmc::update_nconf_counter(hparams.conf_dir, g_heat, i, path_i);
      }
    }
  }

  if (hparams.do_hmc) {
    std::cout << "## Acceptance rate: " << rate / static_cast<double>(hparams.n_meas)
              << std::endl;
    std::string path_final = conf_path_basename + ".final";
    U.save(path_final);
  }

  return (0);
}

void initialize_U_hmc(gaugeconfig<_u1> &U,
                      bool &g_heat,
                      size_t &g_icounter,
                      double &normalisation,
                      const gp::physics &pparams,
                      const gp::hmc_u1 &hparams) {
  if (hparams.restart) {
    std::cout << "## restart " << hparams.restart << "\n";
    std::vector<std::string> v_ncc = io::hmc::read_nconf_counter(hparams.conf_dir);
    g_heat = boost::lexical_cast<bool>(v_ncc[0]);
    g_icounter = std::stoi(v_ncc[1]);
    std::string config_path = v_ncc[2];

    const size_t err = U.load(config_path);
    if (err != 0) {
      std::cout
        << "Error: failed to load initial gauge configuration for hmc. Aborting.\n";
      std::abort();
    }
  } else {
    std::cout << "hotstart " << hparams.seed << " " << hparams.heat << "\n";
    g_heat = (hparams.heat == true) ? 1.0 : 0.0;
    g_icounter = 0;
    hotstart(U, hparams.seed, g_heat);
  }

  double plaquette = flat_spacetime::gauge_energy(U);
  double fac = 2. / U.getndims() / (U.getndims() - 1);
  normalisation = fac / U.getVolume() / double(U.getNc());
  std::cout << "## Initial Plaquette: " << plaquette * normalisation << std::endl;

  random_gauge_trafo(U, 654321);
  plaquette = flat_spacetime::gauge_energy(U);
  std::cout << "## Plaquette after rnd trafo: " << plaquette * normalisation << std::endl;
}
