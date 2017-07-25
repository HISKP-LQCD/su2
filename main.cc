#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"parse_commandline.hh"

#include<iostream>
#include<iomanip>
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
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017)" << endl << endl;  

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

  gaugeconfig U(gparams.Ls, gparams.Lt, gparams.beta);
  if(gparams.restart) {
    err = U.load(gparams.configfilename);
    if(err != 0) {
      return err;
    }
  }
  else {
    U = hotstart(gparams.Ls, gparams.Lt, gparams.seed, gparams.heat);
  }

  double plaquette = gauge_energy(U);
  cout << "Initital Plaquette: " << plaquette/U.getVolume()/N_c/6. << endl; 

  random_gauge_trafo(U, 654321);
  plaquette = gauge_energy(U);
  cout << "Plaquette after rnd trafo: " << plaquette/U.getVolume()/N_c/6. << endl; 

  double rate = 0.;
  for(size_t i = 0; i < gparams.N_meas; i++) {
    std::mt19937 engine(gparams.seed+i);
    rate += sweep(U, engine, delta, N_hit, gparams.beta);
    cout << i << " " << std::scientific << std::setw(15) << gauge_energy(U)/U.getVolume()/N_c/6. << endl;
    if(i > 0 && (i % gparams.N_save) == 0) {
      std::ostringstream os;
      os << "config." << gparams.Ls << "." << gparams.Lt << ".b" << gparams.beta << "." << i << std::ends;
      U.save(os.str());
    }
  }
  cout << rate/static_cast<double>(gparams.N_meas) << endl;

  std::ostringstream os;
  os << "config." << gparams.Ls << "." << gparams.Lt << ".b" << U.getBeta() << ".final" << std::ends;
  U.save(os.str());

  return(0);
}

