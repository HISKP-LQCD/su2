#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"wilsonloop.hh"
#include"md_update.hh"
#include"monomial.hh"
#include"gradient_flow.hh"
#include"energy_density.hh"
#include"parse_commandline.hh"
#include"version.hh"

#include<iostream>
#include<iomanip>
#include<sstream>
#include<vector>
#include<random>
#include<boost/program_options.hpp>

namespace po = boost::program_options;

using std::cout;
using std::endl;


int main(int ac, char* av[]) {
  general_params gparams;
  size_t nstep;
  bool Wloop;
  bool gradient;
  double tmax;

  cout << "## Measuring Tool for SU(2) gauge theory" << endl;
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017)" << endl;
  cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << endl << endl;  

  po::options_description desc("Allowed options");
  add_general_options(desc, gparams);
  // add measure specific options
  desc.add_options()
    ("Wloops", po::value<bool>(&Wloop)->default_value(false), "measure Wilson loops")
    ("gradient", po::value<bool>(&gradient)->default_value(false), "meausre Grandient flow")
    ("nstep", po::value<size_t>(&nstep)->default_value(1), "measure each nstep config")
    ("tmax", po::value<double>(&tmax)->default_value(9.99), "tmax for gradient flow")
    ;

  int err = parse_commandline(ac, av, desc, gparams);
  if(err > 0) {
    return err;
  }

  gaugeconfig U(gparams.Ls, gparams.Lt, gparams.beta);
  //U = hotstart(gparams.Ls, gparams.Lt, 123456, 0.10);

  for(size_t i = gparams.icounter; i < gparams.N_meas*nstep+gparams.icounter; i+=nstep) {
    std::ostringstream os;
    os << "config." << gparams.Ls << "." << gparams.Lt << ".b" << U.getBeta() << "." << i << std::ends;
    U.load(os.str());
    
    double plaquette = gauge_energy(U);
    cout << "## Initital Plaquette: " << plaquette/U.getVolume()/N_c/6. << endl; 
    cout << "## Initial Energy density: " << energy_density(U) << endl;
    
    random_gauge_trafo(U, gparams.seed);
    plaquette = gauge_energy(U);
    cout << "## Plaquette after rnd trafo: " << std::scientific << std::setw(15) << plaquette/U.getVolume()/N_c/6. << endl; 
    cout << "## Energy density: " << energy_density(U) << endl;
    
    if(Wloop) {
      std::ostringstream os;
      os << "wilsonloop.";
      auto prevw = os.width(6);
      auto prevf = os.fill('0');
      os << i;
      os.width(prevw);
      os.fill(prevf);
      os << ".dat" << std::ends;
      compute_special_loops(U, os.str());
    }
    if(gradient) {
      std::ostringstream os;
      os << "gradient_flow.";
      auto prevw = os.width(6);
      auto prevf = os.fill('0');
      os << i;
      os.width(prevw);
      os.fill(prevf);
      gradient_flow(U, os.str(), tmax);
    }
  }

  return(0);
}
