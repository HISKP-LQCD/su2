#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"md_update.hh"
#include"monomial.hh"
#include"integrator.hh"
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
  size_t N_rev;
  size_t n_steps;
  size_t exponent;
  double tau;
  size_t integs;

  cout << "## HMC Algorithm for SU(2) gauge theory" << endl;
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017, 2021)" << endl;
  cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << endl << endl;  

  po::options_description desc("Allowed options");
  add_general_options(desc, gparams);
  // add HMC specific options
  desc.add_options()
    ("nrev", po::value<size_t>(&N_rev)->default_value(0), "frequenz of reversibility tests N_rev, 0: not reversibility test")
    ("nsteps", po::value<size_t>(&n_steps)->default_value(1000), "n_steps")
    ("tau", po::value<double>(&tau)->default_value(1.), "trajectory length tau")
    ("exponent", po::value<size_t>(&exponent)->default_value(0), "exponent for rounding")
    ("integrator", po::value<size_t>(&integs)->default_value(0), "itegration scheme to be used: 0=leapfrog, 1=lp_leapfrog, 2=omf4, 3=lp_omf4")
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
  cout << "## Initital Plaquette: " << plaquette/U.getVolume()/double(U.getNc())/6. << endl; 

  random_gauge_trafo(U, 654321);
  plaquette = gauge_energy(U);
  cout << "## Plaquette after rnd trafo: " << plaquette/U.getVolume()/double(U.getNc())/6. << endl; 

  // Molecular Dynamics parameters
  md_params mdparams(n_steps, tau);
  
  // generate list of monomials
  gaugemonomial<double, su2> gm(0);
  kineticmonomial<double, su2> km(0);
  km.setmdpassive();

  std::list<monomial<double, su2>*> monomial_list;
  monomial_list.push_back(&gm);
  monomial_list.push_back(&km);

  integrator<double, su2> * md_integ = set_integrator<double, su2>(integs, exponent);

  std::ofstream os;
  if(gparams.icounter == 0) 
    os.open("output.hmc.data", std::ios::out);
  else 
    os.open("output.hmc.data", std::ios::app);

  double rate = 0.;
  for(size_t i = gparams.icounter; i < gparams.N_meas+gparams.icounter; i++) {
    mdparams.disablerevtest();
    if(i > 0 && N_rev != 0 && (i) % N_rev == 0) {
      mdparams.enablerevtest();
    }
    // PRNG engine
    std::mt19937 engine(gparams.seed+i);
    // perform the MD update
    md_update(U, engine, mdparams, monomial_list, *md_integ);

    double energy = gauge_energy(U);
    rate += mdparams.getaccept();
    cout << i << " " << mdparams.getaccept() << " " << std::scientific << std::setw(18) << std::setprecision(15) << energy/U.getVolume()/double(U.getNc())/6. << " " << std::setw(15) << mdparams.getdeltaH() << " " 
         << std::setw(15) << rate/static_cast<double>(i+1) << " ";
    if(mdparams.getrevtest()) {
      cout << mdparams.getdeltadeltaH();
    }
    else cout << "NA";
    cout << endl;

    os << i << " " << mdparams.getaccept() << " " << std::scientific << std::setw(18) << std::setprecision(15) << energy/U.getVolume()/double(U.getNc())/6. << " " << std::setw(15) << mdparams.getdeltaH() << " " 
       << std::setw(15) << rate/static_cast<double>(i+1) << " ";
    if(mdparams.getrevtest()) {
      os << mdparams.getdeltadeltaH();
    }
    else os << "NA";
    os << endl;

    if(i > 0 && (i % gparams.N_save) == 0) {
      std::ostringstream oss;
      oss << "config." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "." << gparams.Lt << ".b" << gparams.beta << "." << i << std::ends;
      U.save(oss.str());
    }
  }
  cout << "## Acceptance rate: " << rate/static_cast<double>(gparams.N_meas) << endl;

  std::ostringstream oss;
  oss << "config." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "." << gparams.Lt << ".b" << U.getBeta() << ".final" << std::ends;
  U.save(oss.str());
  return(0);
}
