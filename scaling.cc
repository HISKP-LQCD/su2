#include"su2.hh"
#include"u1.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"parse_commandline.hh"
#include"energy_density.hh"
#include"version.hh"
#include"wilsonloop.hh"

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
#include<chrono>

using std::vector;
using std::cout;
using std::endl;
namespace po = boost::program_options;

int main(int ac, char* av[]) {
  general_params gparams;

  size_t N_hit = 10;
  double delta = 0.1;

  cout << "## Measuring the scaling of parallelization of the U(1) functions" << endl;
  //~ cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017, 2021)" << endl;
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

  gaugeconfig<_u1> U(gparams.Lx, gparams.Ly, gparams.Lz, gparams.Lt, gparams.ndims, gparams.beta);
  
 
  
  double fac = 2./U.getndims()/(U.getndims()-1);
  const double normalisation = fac/U.getVolume();
  size_t facnorm=gparams.ndims>2?gparams.ndims/(gparams.ndims-2):0;

  //set things up for parallel computing in sweep
  #ifdef _USE_OMP_
  int threads=omp_get_max_threads();
  #else
  int threads=1;
  #endif   
  
  char filename[100];
  sprintf(filename, "resultscalingNt%luNs%lubeta%fxi%fmaxthreads%dnmeas%lunsave%lu", gparams.Lt, gparams.Lx, gparams.beta, gparams.xi, threads, gparams.N_meas, gparams.N_save);
  std::ofstream os;
  os.open(filename, std::ios::out);
  os << std::setw(14) << "##threads  " << "time_sweep  " << "speedup_sweep  " << "time_loops  " << "speedup_loops  " << std::endl;
  double rate = 0.;
  std::chrono::duration<double, std::micro> elapse_sweep_one, elapse_loop_one;
  for(size_t measurement=1;measurement<10;measurement++){
  for(size_t thread=1;thread<=threads;thread++){
    hotstart(U, gparams.seed, gparams.heat);
    omp_set_num_threads(thread);
    auto start = std::chrono::high_resolution_clock::now();
    for(size_t i = gparams.icounter; i < gparams.N_meas*thread + gparams.icounter; i+=thread) {
      std::mt19937 * engines =new std::mt19937[thread];
      for(size_t engine=0;engine<thread;engine+=1){
        engines[engine].seed(gparams.seed+i+engine);
      }
      size_t inew = (i-gparams.icounter)/thread+gparams.icounter;//counts loops, loop-variable needed too have one RNG per thread with different seeds 
      rate += sweep(U, engines, delta, N_hit, gparams.beta, gparams.xi, gparams.anisotropic);
      double energy = gauge_energy(U, true);
      double E = 0., Q = 0.;
      energy_density(U, E, Q);
      //measuring spatial plaquettes only means only half of all plaquettes are measured, so need factor 2 for normalization to 1
      cout << inew << " " << std::scientific << std::setw(18) << std::setprecision(15) << energy*normalisation*facnorm << "  " << Q << endl;
      if(inew > 0 && (inew % gparams.N_save) == 0) {
        std::ostringstream oss;
        oss << "configu1." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz<< "." << gparams.Lt << ".b" << std::fixed << gparams.beta << ".x" << gparams.xi << "." << inew << std::ends;
        U.save(oss.str());
      }
      delete engines;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> elapsed_time= end-start;
    if(thread==1){
        elapse_sweep_one=elapsed_time;
    }
    os << thread << "  " << std::setw(14) << std::scientific << elapsed_time.count() << "  " << elapse_sweep_one.count()/elapsed_time.count();
    
      //~ open file for saving results
  std::ofstream resultfile;
  char filenamepot[200];
  
  
  double loop;
  
  start = std::chrono::high_resolution_clock::now();
  for(size_t i = gparams.icounter+gparams.N_save; i < gparams.N_meas+gparams.icounter; i+=gparams.N_save) {
    std::ostringstream oss; 
    oss << "configu1." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "." << gparams.Lt << ".b" << std::fixed << U.getBeta() << ".x" << gparams.xi << "." << i << std::ends;
    U.load(oss.str());
    //~ //calculate wilsonloops
    if(gparams.ndims==4){
    for (size_t x=1 ; x<=gparams.Lx/2 ; x++){
      sprintf(filenamepot, "result.u1potential.Nt%lu.Ns%lu.b%f.xi%f.x%lu",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, x);
      resultfile.open(filenamepot, std::ios::app);
    
    //Measure for two radii each time by changing one of the coordinates not needed for the measurement
    //Measure (x,t) and (x,y), with "t" the anisotropic direction, "x" the "first" isotropic direction and "y" taken as the average of the other two directions
      for (size_t t=1 ; t<=gparams.Lt/2 ; t++){
        loop=wilsonloop_non_planar(U, {t, x, 0, 0});
        loop+=wilsonloop_non_planar(U, {t, 0, x, 0});
        loop+=wilsonloop_non_planar(U, {t, 0, 0, x});
        resultfile << std::setw(14) << std::scientific << loop/U.getVolume()/3.0 << "  " ; 
        
        loop= wilsonloop_non_planar(U, {t, x, 1, 0});
        loop+=wilsonloop_non_planar(U, {t, x, 0, 1});
        loop+=wilsonloop_non_planar(U, {t, 1, x, 0});
        loop+=wilsonloop_non_planar(U, {t, 0, x, 1});
        loop+=wilsonloop_non_planar(U, {t, 0, 1, x});
        loop+=wilsonloop_non_planar(U, {t, 1, 0, x});
        resultfile << std::setw(14) << std::scientific << loop/U.getVolume()/6.0 << "  " ;
      }
      for (size_t y=1 ; y<=gparams.Lx/2 ; y++){ 
        loop= wilsonloop_non_planar(U, {0, x, y, 0});
        loop+=wilsonloop_non_planar(U, {0, x, 0, y});
        loop+=wilsonloop_non_planar(U, {0, y, x, 0});
        loop+=wilsonloop_non_planar(U, {0, 0, x, y});
        loop+=wilsonloop_non_planar(U, {0, 0, y, x});
        loop+=wilsonloop_non_planar(U, {0, y, 0, x});
        resultfile << std::setw(14) << std::scientific << loop/U.getVolume()/6.0 << "  " ;
        
        loop= wilsonloop_non_planar(U, {1, x, y, 0});
        loop+=wilsonloop_non_planar(U, {1, x, 0, y});
        loop+=wilsonloop_non_planar(U, {1, y, x, 0});
        loop+=wilsonloop_non_planar(U, {1, 0, x, y});
        loop+=wilsonloop_non_planar(U, {1, 0, y, x});
        loop+=wilsonloop_non_planar(U, {1, y, 0, x});
        resultfile << std::setw(14) << std::scientific << loop/U.getVolume()/6.0 << "  " ;
      }
    resultfile << std::endl;
    resultfile.close(); 
    }
    }
    }
    
    end = std::chrono::high_resolution_clock::now();
    elapsed_time= end-start;
    if(thread==1){
        elapse_loop_one=elapsed_time;
    }
    os << "  " << std::setw(14) << std::scientific << elapsed_time.count() << "  " << elapse_loop_one.count()/elapsed_time.count();
    os << std::endl;
}
}
 
  os.close();
  return(0);
}

/**
 * set up parallel, resultfile
 * for each number of threads
 * measure nmeas configurations
 * for each configuration measure gauge energy and wilson loops
 * resultfile:
 * nthreads, time(sweeps), time(energy), time(wilson)
 * **/

