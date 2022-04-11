#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"parse_commandline.hh"
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

using std::vector;
using std::cout;
using std::endl;
namespace po = boost::program_options;

int main(int ac, char* av[]) {
  general_params gparams;

  size_t N_hit = 10;
  double delta = 0.1;

  cout << "## Metropolis Algorithm for SU(2) gauge theory" << endl;
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017,2020)" << endl;
  cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << endl << endl;  

  po::options_description desc("Allowed options");
  add_general_options(desc, gparams);

  // add Metropolis specific options
  desc.add_options()
    ("nhit", po::value<size_t>(&N_hit)->default_value(10), "N_hit")
    ("delta,d", po::value<double>(&delta), "delta")
    ;

  int err = parse_commandline(ac, av, desc, gparams);
  if(err > 0) {
    return err;
  }
  
  #ifdef _USE_OMP_
  bool parallel = true;
  #else
  bool parallel = false;
  #endif
  /**
   * the parallelisation of the sweep-function first iterates over all odd points in t and then over all even points
   * because the nearest neighbours must not change during the updates, this is not possible for an uneven number of points in T
   * */
  if (gparams.Lt%2 != 0 && parallel){
    std::cerr << "For parallel computing an even number of points in T is needed!" << std::endl;
    omp_set_num_threads(1);
    std::cerr << "Continuing with one thread." << std::endl;
  }    
  
  // load/set initial configuration    
  gaugeconfig<su2> U(gparams.Lx, gparams.Ly, gparams.Lz, gparams.Lt, gparams.ndims, gparams.beta);
  if(gparams.restart) {
    err = U.load(gparams.configfilename);
    if(err != 0) {
      return err;
    }
  }
  else {
    hotstart(U, gparams.seed, gparams.heat);
  }
  
  // check gauge invariance, set up factors needed to normalise plaquette, spacial plaquette
  double plaquette = gauge_energy(U);
  double fac = 2./U.getndims()/(U.getndims()-1);
  const double normalisation = fac/U.getVolume()/double(U.getNc());
  cout << "Initital Plaquette: " << plaquette*normalisation << endl; 

  random_gauge_trafo(U, 654321);
  plaquette = gauge_energy(U);
  cout << "Plaquette after rnd trafo: " << plaquette*normalisation << endl; 

  //set things up for parallel computing in sweep
  #ifdef _USE_OMP_
  int threads=omp_get_max_threads();
  #else
  int threads=1;
  #endif    

  std::ofstream os;
  std::ofstream acceptancerates;
  if(gparams.icounter == 0) 
    os.open("output.metropolis.data", std::ios::out);
  else
    os.open("output.metropolis.data", std::ios::app);
  std::vector<double> rate = {0., 0.};
  
  /**
   * do measurements:
   * sweep: do N_hit Metropolis-Updates of every link in the lattice
   * calculate plaquette, spacial plaquette, energy density and write to stdout and output-file
   * save every nave configuration
   * */  
  for(size_t i = gparams.icounter; i < gparams.n_meas*threads + gparams.icounter; i+=threads) {
    std::vector<std::mt19937> engines(threads);
    for(size_t engine = 0 ; engine < threads ; engine += 1){
      engines[engine].seed(gparams.seed + i + engine);
    }
    //inew counts loops, loop-variable needed to have one RNG per thread with different seeds for every measurement
    size_t inew = (i-gparams.icounter) / threads + gparams.icounter;
    rate += sweep(U, engines, delta, N_hit, gparams.beta);
    double energy = gauge_energy(U);
    double E = 0., Q = 0.;
    energy_density(U, E, Q);
    cout << inew << " " << std::scientific << std::setw(18) << std::setprecision(15) << energy*normalisation << " " << Q << endl;
    os << inew << " " << std::scientific << std::setw(18) << std::setprecision(15) << energy*normalisation << " " << Q << endl;
    if(inew > 0 && (inew % gparams.N_save) == 0) {
      std::ostringstream oss;
      oss << "config." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "." << gparams.Lt << ".b" << gparams.beta << "." << inew << std::ends;
      U.save(oss.str());
    }
  }
  // save acceptance rates to additional file to keep track of measurements
  cout << "## Acceptance rate " << rate[0]/static_cast<double>(gparams.n_meas) << " temporal acceptance rate " << rate[1]/static_cast<double>(gparams.n_meas) << endl;
  acceptancerates.open("acceptancerates.data", std::ios::app);
  acceptancerates << rate[0]/static_cast<double>(gparams.n_meas) << " " << rate[1]/static_cast<double>(gparams.n_meas) << " "
   << gparams.beta << " " << gparams.Lx << " " << gparams.Lt << " " << gparams.xi << " " 
   << delta << " " << gparams.heat << " " << threads << " " << N_hit << " " << gparams.n_meas << " " << gparams.seed << " " << endl;
  acceptancerates.close();

  std::ostringstream oss;
  oss << "config." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "." << gparams.Lt << ".b" << U.getBeta() << ".final" << std::ends;
  U.save(oss.str());

  return(0);
}

