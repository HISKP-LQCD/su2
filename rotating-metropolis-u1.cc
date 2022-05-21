/**
 * @file rotating-metropolis-u1.cc
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief
 * @version 0.1
 * @date 2022-05-20
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "rotating-energy_density.hpp"
#include "rotating-gauge_energy.hpp"
#include "gaugeconfig.hh"
#include "output.hh"
#include "parse_input_file.hh"
#include "random_gauge_trafo.hh"
#include "su2.hh"
#include "sweep.hh"
#include "u1.hh"
#include "vectorfunctions.hh"
#include "version.hh"

#ifdef _USE_OMP_
#include <omp.h>
#endif

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

namespace po = boost::program_options;

int main(int ac, char *av[]) {
  std::cout << "## Metropolis Algorithm for U(1) gauge theory" << std::endl;
  std::cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH
            << std::endl
            << std::endl;


  namespace gp = global_parameters;
  gp::physics pparams; // physics parameters
  gp::metropolis_u1 mcparams; // mcmc parameters

  std::string input_file; // yaml input file path
  int err = input_file_parsing::parse_command_line(ac, av, input_file);
  if (err > 0) {
    return err;
  }

  namespace in_hmc = input_file_parsing::u1::hmc;

  err = in_hmc::parse_input_file(input_file, pparams, mcparams);
  if (err > 0) {
    return err;
  }

  boost::filesystem::create_directories(boost::filesystem::absolute(mcparams.conf_dir));

#ifdef _USE_OMP_
  /**
   * the parallelisation of the sweep-function first iterates over all odd points in t and
   * then over all even points because the nearest neighbours must not change during the
   * updates, this is not possible for an uneven number of points in T
   * */
  if (pparams.Lt % 2 != 0) {
    std::cerr << "For parallel computing an even number of points in T is needed!"
              << std::endl;
    omp_set_num_threads(1);
    std::cerr << "Continuing with one thread." << std::endl;
  }
  // set things up for parallel computing in sweep
  int threads = omp_get_max_threads();
#else
  int threads = 1;
#endif
  // std::cout << "threads " << threads << std::endl;

  // load/set initial configuration
  gaugeconfig<_u1> U(pparams.Lx, pparams.Ly, pparams.Lz, pparams.Lt, pparams.ndims,
                     pparams.beta);
  if (mcparams.restart) {
    std::cout << "restart " << mcparams.restart << std::endl;
    err = U.load(mcparams.configfilename);
    if (err != 0) {
      return err;
    }
  } else {
    hotstart(U, mcparams.seed, mcparams.heat);
  }

  // check gauge invariance, set up factors needed to normalise plaquette, spatial
  // plaquette
  double plaquette = rotating_spacetime::gauge_energy(U, pparams.Omega);
  double fac = 2. / U.getndims() / (U.getndims() - 1);
  const double normalisation = fac / U.getVolume();
  size_t facnorm = (pparams.ndims > 2) ? pparams.ndims / (pparams.ndims - 2) : 0;

  std::cout << "## Initial Plaquette: " << plaquette * normalisation << std::endl;

  random_gauge_trafo(U, 654321);
  plaquette = rotating_spacetime::gauge_energy(U, pparams.Omega);
  std::cout << "## Plaquette after rnd trafo: " << plaquette * normalisation << std::endl;

  std::ofstream os;
  std::ofstream acceptancerates;
  if (mcparams.icounter == 0)
    os.open(mcparams.conf_dir + "/output.u1-metropolis.data", std::ios::out);
  else
    os.open(mcparams.conf_dir + "/output.u1-metropolis.data", std::ios::app);
  std::vector<double> rate = {0., 0.};

  // get basename for configs
  std::string conf_path_basename = output::get_conf_path_basename(pparams, mcparams);

  /**
   * do measurements:
   * sweep: do N_hit Metropolis-Updates of every link in the lattice
   * calculate plaquette, spacial plaquette, energy density with and without cloverdef and
   * write to stdout and output-file save every nave configuration
   * */
  for (size_t i = mcparams.icounter; i < mcparams.n_meas * threads + mcparams.icounter;
       i += threads) {
    std::vector<std::mt19937> engines(threads);
    for (size_t engine = 0; engine < threads; engine += 1) {
      engines[engine].seed(mcparams.seed + i + engine);
    }
    // inew counts loops, loop-variable needed to have one RNG per thread with different
    // seeds for every measurement
    size_t inew = (i - mcparams.icounter) / threads + mcparams.icounter;
    rate += sweep(U, engines, mcparams.delta, mcparams.N_hit, pparams.beta, pparams.xi,
                  pparams.anisotropic);

    double energy = rotating_spacetime::gauge_energy(U, pparams.Omega, true);
    double E = 0., Q = 0.;
    rotating_spacetime::energy_density(U, pparams.Omega, E, Q);
    // measuring spatial plaquettes only means only (ndims-1)/ndims of all plaquettes are
    // measured, so need facnorm for normalization to 1
    std::cout << inew << " " << std::scientific << std::setw(18) << std::setprecision(15)
              << energy * normalisation * facnorm << " ";
    os << inew << " " << std::scientific << std::setw(18) << std::setprecision(15)
       << energy * normalisation * facnorm << " ";

    energy = rotating_spacetime::gauge_energy(U, pparams.Omega, false);
    std::cout << energy * normalisation << " " << Q << " ";
    os << energy * normalisation << " " << Q << " ";
    rotating_spacetime::energy_density(U, pparams.Omega, E, Q, false);
    std::cout << Q << std::endl;
    os << Q << std::endl;

    if (inew > 0 && (inew % mcparams.N_save) == 0) {
      std::ostringstream oss_i;
      oss_i << conf_path_basename << "." << inew << std::ends;
      U.save(mcparams.conf_dir + "/" + oss_i.str());
    }
  }
  // save acceptance rates to additional file to keep track of measurements
  std::cout << "## Acceptance rate " << rate[0] / static_cast<double>(mcparams.n_meas)
            << " temporal acceptance rate "
            << rate[1] / static_cast<double>(mcparams.n_meas) << std::endl;
  acceptancerates.open(mcparams.conf_dir + "/acceptancerates.data", std::ios::app);
  acceptancerates << rate[0] / static_cast<double>(mcparams.n_meas) << " "
                  << rate[1] / static_cast<double>(mcparams.n_meas) << " " << pparams.beta
                  << " " << pparams.Lx << " " << pparams.Lt << " " << pparams.xi << " "
                  << mcparams.delta << " " << mcparams.heat << " " << threads << " "
                  << mcparams.N_hit << " " << mcparams.n_meas << " " << mcparams.seed
                  << " " << std::endl;
  acceptancerates.close();

  std::ostringstream oss;
  oss << conf_path_basename << ".final" << std::ends;
  U.save(mcparams.conf_dir + "/" + oss.str());
  return (0);
}
