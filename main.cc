#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"parse_commandline.hh"
#include"version.hh"

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
  double fac = 1.;
  if(U.getndims() == 4) fac = 1./6.;
  if(U.getndims() == 3) fac = 1./2.;
  const double normalisation = fac/U.getVolume()/N_c;
  cout << "Initital Plaquette: " << plaquette*normalisation << endl; 

  random_gauge_trafo(U, 654321);
  plaquette = gauge_energy(U);
  cout << "Plaquette after rnd trafo: " << plaquette*normalisation << endl; 

  std::ofstream os;
  if(gparams.icounter == 0) 
    os.open("output.metropolis.data", std::ios::out);
  else
    os.open("output.metropolis.data", std::ios::app);
  double rate = 0.;
  for(size_t i = gparams.icounter; i < gparams.N_meas + gparams.icounter; i++) {
    std::mt19937 engine(gparams.seed+i);
    rate += sweep(U, engine, delta, N_hit, gparams.beta);
    double energy = gauge_energy(U);
    cout << i << " " << std::scientific << std::setw(18) << std::setprecision(15) << energy*normalisation << " " << -U.getBeta()/N_c*(U.getVolume()*N_c/fac - energy) << endl;
    os << i << " " << std::scientific << std::setw(18) << std::setprecision(15) << energy*normalisation << " " << -U.getBeta()/N_c*(U.getVolume()*N_c/fac - energy) << endl;
    if(i > 0 && (i % gparams.N_save) == 0) {
      std::ostringstream oss;
      oss << "config." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "." << gparams.Lt << ".b" << gparams.beta << "." << i << std::ends;
      U.save(oss.str());
    }
  }
  cout << rate/static_cast<double>(gparams.N_meas) << endl;

  std::ostringstream oss;
  oss << "config." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "." << gparams.Lt << ".b" << U.getBeta() << ".final" << std::ends;
  U.save(oss.str());

  return(0);
}

