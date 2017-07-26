#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"kramers_md_update.hh"
#include"monomial.hh"
#include"parse_commandline.hh"
#include"version.hh"

#include<iostream>
#include<fstream>
#include<iomanip>
#include<sstream>
#include<random>
#include<boost/program_options.hpp>

using std::cout;
using std::endl;
namespace po = boost::program_options;

int main(int ac, char* av[]) {
  const size_t n_steps = 1;

  general_params gparams;
  size_t N_rev;
  size_t exponent;
  double tau;

  cout << "## Kramers Algorithm for SU(2) gauge theory" << endl;
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017)" << endl;
  cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << endl << endl;  

  po::options_description desc("Allowed options");
  add_general_options(desc, gparams);
  // add HMC specific options
  desc.add_options()
    ("tau", po::value<double>(&tau)->default_value(1.), "trajectory length tau")
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

  mdparams.setkmax(5);

  std::ofstream os;
  if(gparams.icounter == 0) 
    os.open("output.kramers.data", std::ios::out);
  else 
    os.open("output.kramers.data", std::ios::app);

  double rate = 0.;
  for(size_t i = gparams.icounter; i < gparams.N_meas + gparams.icounter; i++) {
    mdparams.disablerevtest();
    kramers_md_update(U, engine, mdparams, monomial_list);

    double energy = gauge_energy(U);
    rate += mdparams.getaccept();
    cout << i << " " << mdparams.getaccept() << " " << std::scientific << std::setw(18)  << std::setprecision(15) << energy/U.getVolume()/N_c/6. << " " << std::setw(15) << mdparams.getdeltaH() << " " 
         << std::setw(15) << rate/static_cast<double>(i+1) << std::endl;

    os << i << " " << mdparams.getaccept() << " " << std::scientific << std::setw(18)  << std::setprecision(15) << energy/U.getVolume()/N_c/6. << " " << std::setw(15) << mdparams.getdeltaH() << " " 
       << std::setw(15) << rate/static_cast<double>(i+1) << std::endl;

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
