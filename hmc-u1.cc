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

#include "energy_density.hh"
#include "flat-gauge_energy.hpp"
#include "gaugeconfig.hh"
#include "integrator.hh"
#include "md_update.hh"
#include "monomial.hh"
#include "omeasurements.hpp"
#include "output.hh"
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

int main(int ac, char *av[]) {
  std::cout << "## HMC Algorithm for U(1) gauge theory" << std::endl;
  std::cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017, 2021)"
            << std::endl;
  std::cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH
            << std::endl;

  namespace gp = global_parameters;
  gp::physics pparams; // physics parameters
  gp::hmc_u1 hparams; // hmc parameters

  std::string input_file; // yaml input file path
  int err = input_file_parsing::parse_command_line(ac, av, input_file);
  if (err > 0) {
    return err;
  }

  namespace in_hmc = input_file_parsing::u1::hmc;

  err = in_hmc::parse_input_file(input_file, pparams, hparams);
  if (err > 0) {
    return err;
  }

  boost::filesystem::create_directories(boost::filesystem::absolute(hparams.conf_dir));

  gaugeconfig<_u1> U(pparams.Lx, pparams.Ly, pparams.Lz, pparams.Lt, pparams.ndims,
                     pparams.beta);
  if (hparams.restart) {
    std::cout << "restart " << hparams.restart << "\n";
    err = U.load(hparams.configfilename);
    if (err != 0) {
      return err;
    }
  } else {
    std::cout << "hotstart " << hparams.seed << " " << hparams.heat << "\n";
    const double heat_val = (hparams.heat == true) ? 1.0 : 0.0;
    hotstart(U, hparams.seed, heat_val);
  }

  double plaquette = flat_spacetime::gauge_energy(U);
  double fac = 2. / U.getndims() / (U.getndims() - 1);
  const double normalisation = fac / U.getVolume() / double(U.getNc());
  std::cout << "## Initital Plaquette: " << plaquette * normalisation << std::endl;

  random_gauge_trafo(U, 654321);
  plaquette = flat_spacetime::gauge_energy(U);
  std::cout << "## Plaquette after rnd trafo: " << plaquette * normalisation << std::endl;

  // Molecular Dynamics parameters
  md_params mdparams(hparams.n_steps, hparams.tau);

  // generate list of monomials
  std::list<monomial<double, _u1> *> monomial_list;
  gaugemonomial<double, _u1> gm(0);
  rotating_spacetime::gauge_monomial<double, _u1> gm_rot(0, pparams.Omega);

  kineticmonomial<double, _u1> km(0);
  km.setmdpassive();
  monomial_list.push_back(&km);

  staggered::detDDdag_monomial<double, _u1> detDDdag(
    0, pparams.m0, hparams.solver, hparams.tolerance_cg, hparams.seed_pf,
    hparams.solver_verbosity);

  monomial_list.push_back(&gm);

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
  if (hparams.icounter == 0)
    os.open(hparams.conf_dir + "/output.hmc.data", std::ios::out);
  else
    os.open(hparams.conf_dir + "/output.hmc.data", std::ios::app);

  std::cout << "## Normalization factor: A = 2/(d*(d-1)*N_lat*N_c) = " << std::scientific
            << std::setw(18) << std::setprecision(15) << normalisation << "\n";
  std::cout << "## Acceptance rate parcentage: rho = rate/(i+1)\n";

  // header: column names in the output
  std::string head_str = output::hmc::get_header(" ");
  std::cout << head_str;
  os << head_str;

  double rate = 0.;

  const std::string conf_path_basename = output::get_conf_path_basename(pparams, hparams);

  for (size_t i = hparams.icounter; i < hparams.n_meas + hparams.icounter; i++) {
    mdparams.disablerevtest();
    if (i > 0 && hparams.N_rev != 0 && (i) % hparams.N_rev == 0) {
      mdparams.enablerevtest();
    }
    // PRNG engine
    std::mt19937 engine(hparams.seed + i);
    // perform the MD update
    md_update(U, engine, mdparams, monomial_list, *md_integ);

    double energy = flat_spacetime::gauge_energy(U);
    double E = 0., Q = 0.;
    energy_density(U, E, Q);
    rate += mdparams.getaccept();

    std::cout << i << " " << mdparams.getaccept() << " " << std::scientific
              << std::setw(18) << std::setprecision(15) << energy * normalisation << " "
              << std::setw(15) << mdparams.getdeltaH() << " " << std::setw(15)
              << rate / static_cast<double>(i + 1) << " ";

    if (mdparams.getrevtest()) {
      std::cout << mdparams.getdeltadeltaH();
    } else
      std::cout << "NA";
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

    if (i > 0 && (i % hparams.N_save) == 0) { // saving U after each N_save trajectories
      std::string path_i = conf_path_basename + "." + std::to_string(i);
      U.save(path_i);
      if (hparams.make_omeas && mdparams.getaccept() && i > hparams.omeas.icounter &&
          i % hparams.omeas.nstep == 0) {
        if (hparams.omeas.Wloop) {
          if (hparams.omeas.verbosity > 0) {
            std::cout << "online measuring Wilson loop\n";
          }
          omeasurements::meas_wilson_loop<_u1>(U, i, hparams.conf_dir);
        }
        if (hparams.omeas.gradient) {
          if (hparams.omeas.verbosity > 0) {
            std::cout << "online measuring: Gradient flow";
          }
          omeasurements::meas_gradient_flow<_u1>(U, i, hparams.conf_dir,
                                                 hparams.omeas.tmax);
        }

        if (hparams.omeas.pion_staggered) {
          if (hparams.omeas.verbosity > 0) {
            std::cout << "online measuring Pion correlator\n";
          }
          omeasurements::meas_pion_correlator<_u1>(U, i, pparams.m0, hparams);
        }
      }
    }
  }
  std::cout << "## Acceptance rate: " << rate / static_cast<double>(hparams.n_meas)
            << std::endl;

  std::string path_final = conf_path_basename + ".final";
  U.save(path_final);

  return (0);
}
