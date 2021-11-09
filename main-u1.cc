#include"su2.hh"
#include"u1.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"parse_commandline.hh"
#include"energy_density.hh"
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
  double xi = 1.0;

  cout << "## Metropolis Algorithm for U(1) gauge theory" << endl;
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017, 2021)" << endl;
  cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << endl << endl;  

  po::options_description desc("Allowed options");
  add_general_options(desc, gparams);

  // add Metropolis specific options
  desc.add_options()
    ("nhit", po::value<size_t>(&N_hit)->default_value(10), "N_hit")
    ("delta,d", po::value<double>(&delta), "delta")
    ("tau", po::value<double>(&xi)->default_value(1.0), "xi, characteristic of anisotropy")
    ;

  int err = parse_commandline(ac, av, desc, gparams);
  if(err > 0) {
    return err;
  }

  gaugeconfig<_u1> U(gparams.Lx, gparams.Ly, gparams.Lz, gparams.Lt, gparams.ndims, gparams.beta);
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
  const double normalisation = fac/U.getVolume();
  cout << "## Initital Plaquette: " << plaquette*normalisation << endl; 

  random_gauge_trafo(U, 654321);
  plaquette = gauge_energy(U);
  cout << "## Plaquette after rnd trafo: " << plaquette*normalisation << endl; 


  std::ofstream os;
  if(gparams.icounter == 0) 
    os.open("output.u1-metropolis.data", std::ios::out);
  else
    os.open("output.u1-metropolis.data", std::ios::app);
  double rate = 0.;
  for(size_t i = gparams.icounter; i < gparams.N_meas + gparams.icounter; i++) {
    std::mt19937 engine(gparams.seed+i);
    rate += sweep(U, engine, delta, N_hit, gparams.beta, xi);
    double energy = gauge_energy(U);
    double E = 0., Q = 0.;
    //~ energy_density(U, E, Q);
    cout << i << " " << std::scientific << std::setw(18) << std::setprecision(15) << energy*normalisation << " " << Q << " ";
    os << i << " " << std::scientific << std::setw(18) << std::setprecision(15) << energy*normalisation << " " << Q << " ";
    //~ energy_density(U, E, Q, false);
    cout << Q << endl;
    os << Q << endl;
    if(i > 0 && (i % gparams.N_save) == 0) {
      std::ostringstream oss;
      oss << "configu1." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz<< "." << gparams.Lt << ".b" << gparams.beta << "." << i << std::ends;
      U.save(oss.str());
    }
  }
  cout << "## Acceptance rate " << rate/static_cast<double>(gparams.N_meas) << endl;

  std::ostringstream oss;
  oss << "configu1." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz<< "." << gparams.Lt << ".b" << U.getBeta() << ".final" << std::ends;
  U.save(oss.str());

  return(0);
}

