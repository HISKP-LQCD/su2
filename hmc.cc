#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"md_update.hh"
#include"monomial.hh"
#include"parse_commandline.hh"
#include"version.hh"

#include<iostream>
#include<fstream>
#include<iomanip>
#include<sstream>
#include<random>
#include<boost/program_options.hpp>

namespace po = boost::program_options;
using std::cout;
using std::endl;

int main(int ac, char* av[]) {
  general_params gparams;
  size_t N_rev = 1;
  size_t n_steps = 10;
  size_t exponent = 16;
  double tau = 1.;

  cout << "## HMC Algorithm for SU(2) gauge theory" << endl;
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017)" << endl;
  cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << endl << endl;  

  po::options_description desc("Allowed options");
  add_general_options(desc, gparams);
  // add HMC specific options
  desc.add_options()
    ("nrev", po::value<size_t>(&N_rev)->default_value(0), "frequenz of reversibility tests N_rev, 0: not reversibility test")
    ("nsteps", po::value<size_t>(&n_steps)->default_value(1000), "n_steps")
    ("tau", po::value<double>(&tau)->default_value(1.), "trajectory length tau")
    ("exponent", po::value<size_t>(&exponent)->default_value(0), "exponent for rounding")
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
  // Molecular Dynamics parameters
  md_params mdparams(n_steps, tau);
  // PRNG engine
  std::mt19937 engine(gparams.seed);

  double plaquette = gauge_energy(U);
  cout << "## Initital Plaquette: " << plaquette/U.getVolume()/N_c/6. << endl; 

  random_gauge_trafo(U, 654321);
  plaquette = gauge_energy(U);
  cout << "## Plaquette after rnd trafo: " << plaquette/U.getVolume()/N_c/6. << endl; 

  // generate list of monomials
  gaugemonomial<double> gm(0);
  kineticmonomial<double> km(0);
  km.setmdpassive();

  std::list<monomial<double>*> monomial_list;
  monomial_list.push_back(&gm);
  monomial_list.push_back(&km);

  mdparams.setexponent(exponent);
  std::ofstream os("output.hmc.data", std::ios::app);
  double rate = 0.;
  for(size_t i = gparams.icounter; i < gparams.N_meas+gparams.icounter; i++) {
    mdparams.disablerevtest();
    if(i > 0 && N_rev != 0 && (i) % N_rev == 0) {
      mdparams.enablerevtest();
    }
    md_update(U, engine, mdparams, monomial_list);

    double energy = gauge_energy(U);
    rate += mdparams.getaccept();
    cout << i << " " << mdparams.getaccept() << " " << std::scientific << std::setw(18) << std::setprecision(15) << energy/U.getVolume()/N_c/6. << " " << std::setw(15) << mdparams.getdeltaH() << " " 
         << std::setw(15) << rate/static_cast<double>(i+1) << " ";
    if(mdparams.getrevtest()) {
      cout << mdparams.getdeltadeltaH();
    }
    else cout << "NA";
    cout << endl;

    os << i << " " << mdparams.getaccept() << " " << std::scientific << std::setw(18) << std::setprecision(15) << energy/U.getVolume()/N_c/6. << " " << std::setw(15) << mdparams.getdeltaH() << " " 
       << std::setw(15) << rate/static_cast<double>(i+1) << " ";
    if(mdparams.getrevtest()) {
      os << mdparams.getdeltadeltaH();
    }
    else os << "NA";
    os << endl;

    if(i > 0 && (i % gparams.N_save) == 0) {
      std::ostringstream oss;
      oss << "config." << gparams.Ls << "." << gparams.Lt << ".b" << gparams.beta << "." << i << std::ends;
      U.save(oss.str());
    }
  }
  cout << "## Acceptance rate: " << rate/static_cast<double>(gparams.N_meas) << endl;

  std::ostringstream oss;
  oss << "config." << gparams.Ls << "." << gparams.Lt << ".b" << U.getBeta() << ".final" << std::ends;
  U.save(oss.str());
  return(0);
}
