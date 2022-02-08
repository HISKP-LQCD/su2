#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"parse_commandline.hh"
#include"energy_density.hh"
#include"version.hh"

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
  if (gparams.Lt%2 != 0 && parallel){
    std::cerr << "For parallel computing an even number of points in T is needed!" << std::endl;
    omp_set_num_threads(1);
    std::cerr << "Continuing with one thread." << std::endl;
  }    
      
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
  if(gparams.icounter == 0) 
    os.open("output.metropolis.data", std::ios::out);
  else
    os.open("output.metropolis.data", std::ios::app);
  double rate = 0.;
  for(size_t i = gparams.icounter; i < gparams.N_meas*threads + gparams.icounter; i+=threads) {
    std::vector<std::mt19937> engines(threads);
    for(size_t engine=0;engine<threads;engine+=1){
      engines[engine].seed(gparams.seed+i+engine);
    }
    size_t inew = (i-gparams.icounter)/threads+gparams.icounter;//counts loops, loop-variable needed too have one RNG per thread with different seeds 
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
  cout << "## Acceptance rate " << rate/static_cast<double>(gparams.N_meas) << endl;

  std::ostringstream oss;
  oss << "config." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "." << gparams.Lt << ".b" << U.getBeta() << ".final" << std::ends;
  U.save(oss.str());

  return(0);
}

