// main-u1.cc
/**
 * @file main-u1.cc
 * @author Carsten Urbach (urbach@hiskp.uni-bonn.de)
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief Metropolis Algorithm for U(1) gauge theory
 * @version 0.1
 * @date 2022-05-09
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include"su2.hh"
#include"u1.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"parse_input_file.hh"
#include"energy_density.hh"
#include"version.hh"
#include"vectorfunctions.hh"

#ifdef _USE_OMP_
#  include<omp.h>
#endif

#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<vector>
#include<random>
#include<boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

int main(int ac, char* av[]) {
    
  std::cout << "## Metropolis Algorithm for U(1) gauge theory" << std::endl;
//  std::cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017, 2022)" << std::endl;
  std::cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << std::endl << std::endl;  
  
  namespace gp = global_parameters;
  gp::physics pparams; // physics parameters
  gp::metropolis_u1 mcparams; // mcmc parameters

  std::string input_file; // yaml input file path
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce this help message")(
    "file,f", po::value<std::string>(&input_file)->default_value("NONE"),
    "yaml input file");

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 0;
  }
  
  namespace in_metropolis = input_file_parsing::u1::metropolis;
  int err = in_metropolis::parse_input_file(input_file, pparams, mcparams);
  if (err > 0) {
    return 1;
  }

  boost::filesystem::create_directories(boost::filesystem::absolute(mcparams.outdir));

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
  gaugeconfig<_u1> U(pparams.Lx, pparams.Ly, pparams.Lz, pparams.Lt, pparams.ndims, pparams.beta);
  if(mcparams.restart) {
    std::cout << "restart " << mcparams.restart << std::endl;
    err = U.load(mcparams.configfilename);
    if(err != 0) {
      return err;
    }
  }
  else {
    hotstart(U, mcparams.seed, mcparams.heat);
  }
  
  // check gauge invariance, set up factors needed to normalise plaquette, spacial plaquette
  double plaquette = gauge_energy(U);
  double fac = 2./U.getndims()/(U.getndims()-1);
  const double normalisation = fac/U.getVolume();
  size_t facnorm = (pparams.ndims > 2) ? pparams.ndims/(pparams.ndims-2) : 0;
  
  std::cout << "## Initital Plaquette: " << plaquette*normalisation << std::endl; 

  random_gauge_trafo(U, 654321);
  plaquette = gauge_energy(U);
  std::cout << "## Plaquette after rnd trafo: " << plaquette*normalisation << std::endl; 

  
  std::ofstream os;
  std::ofstream acceptancerates;
  if(mcparams.icounter == 0) 
    os.open(mcparams.outdir+"/output.u1-metropolis.data", std::ios::out);
  else
    os.open(mcparams.outdir+"/output.u1-metropolis.data", std::ios::app);
  std::vector<double> rate = {0., 0.};
  
  //set up name for configs
  const std::string conf_basename = mcparams.conf_basename;
  std::stringstream ss_basename;
  ss_basename << mcparams.conf_basename << ".";
  ss_basename << pparams.Lx << "." << pparams.Ly << "." << pparams.Lz << "."
              << pparams.Lt;
  ss_basename << ".b" << std::fixed << std::setprecision(mcparams.beta_str_width)
              << pparams.beta;
  if(pparams.anisotropic){
    ss_basename << ".x" << std::fixed << std::setprecision(mcparams.beta_str_width)
                << pparams.xi;
  }

  /**
   * do measurements:
   * sweep: do N_hit Metropolis-Updates of every link in the lattice
   * calculate plaquette, spacial plaquette, energy density with and without cloverdef and write to stdout and output-file
   * save every nave configuration
   * */
  for(size_t i = mcparams.icounter; i < mcparams.n_meas*threads + mcparams.icounter; i+=threads) {
    std::vector<std::mt19937> engines(threads);
    for(size_t engine = 0 ; engine < threads ; engine += 1){
      engines[engine].seed(mcparams.seed + i + engine);
    }
    //inew counts loops, loop-variable needed to have one RNG per thread with different seeds for every measurement
    size_t inew = (i - mcparams.icounter) / threads + mcparams.icounter;
    rate += sweep(U, engines, mcparams.delta, mcparams.N_hit, pparams.beta, pparams.xi, pparams.anisotropic);
    
    double energy = gauge_energy(U, true);
    double E = 0., Q = 0.;
    energy_density(U, E, Q);
    //measuring spatial plaquettes only means only (ndims-1)/ndims of all plaquettes are measured, so need facnorm for normalization to 1
    std::cout << inew << " " << std::scientific << std::setw(18) << std::setprecision(15) << energy*normalisation*facnorm << " " ;
    os << inew << " " << std::scientific << std::setw(18) << std::setprecision(15) << energy*normalisation*facnorm << " " ;
    
    energy=gauge_energy(U, false);
    std::cout << energy*normalisation << " " << Q << " ";
    os << energy*normalisation << " " << Q << " ";
    energy_density(U, E, Q, false);
    std::cout << Q << std::endl;
    os << Q << std::endl;
    
    if(inew > 0 && (inew % mcparams.N_save) == 0) {
      std::ostringstream oss_i;
      oss_i << ss_basename.str() << "." << inew << std::ends;
      U.save(mcparams.outdir + "/" + oss_i.str());
    }
  }
  // save acceptance rates to additional file to keep track of measurements
  std::cout << "## Acceptance rate " << rate[0]/static_cast<double>(mcparams.n_meas) 
    << " temporal acceptance rate " << rate[1]/static_cast<double>(mcparams.n_meas) << std::endl;
  acceptancerates.open(mcparams.outdir+"/acceptancerates.data", std::ios::app);
  acceptancerates << rate[0]/static_cast<double>(mcparams.n_meas) << " " << rate[1]/static_cast<double>(mcparams.n_meas) << " "
   << pparams.beta << " " << pparams.Lx << " " << pparams.Lt << " " << pparams.xi << " " 
   << mcparams.delta << " " << mcparams.heat << " " << threads << " "
     << mcparams.N_hit << " " << mcparams.n_meas << " " << mcparams.seed << " " << std::endl;
  acceptancerates.close();

std::ostringstream oss;
  oss << ss_basename.str() << ".final" << std::ends;
  U.save(mcparams.outdir + "/" + oss.str());
  return(0);
}

