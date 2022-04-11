#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"wilsonloop.hh"
#include"md_update.hh"
#include"monomial.hh"
#include"gradient_flow.hh"
#include"energy_density.hh"
#include"parse_commandline.hh"
#include"lyapunov.hh"
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
  bool lyapunov;
  double tmax;
  size_t n_steps;
  size_t exponent;
  double tau;
  size_t integs;

  cout << "## Measuring Tool for SU(2) gauge theory" << endl;
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017)" << endl;
  cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << endl << endl;  

  po::options_description desc("Allowed options");
  add_general_options(desc, gparams);
  // add measure specific options
  desc.add_options()
    ("Wloops", po::value<bool>(&Wloop)->default_value(false), "measure Wilson loops")
    ("gradient", po::value<bool>(&gradient)->default_value(false), "meausre Grandient flow")
    ("lyapunov", po::value<bool>(&lyapunov)->default_value(false), "meausre Lyapunov Exponent")
    ("nstep", po::value<size_t>(&nstep)->default_value(1), "measure each nstep config")
    ("tmax", po::value<double>(&tmax)->default_value(9.99), "tmax for gradient flow")
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

  for(size_t i = gparams.icounter; i < gparams.n_meas*nstep+gparams.icounter; i+=nstep) {
    std::ostringstream os;
    os << "config." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "." << gparams.Lt << ".b" << U.getBeta() << "." << i << std::ends;
    U.load(os.str());
    
    double plaquette = gauge_energy(U);
    double density = 0., Q=0.;
    energy_density(U, density, Q);
    cout << "## Initital Plaquette: " << plaquette/U.getVolume()/double(U.getNc())/6. << endl; 
    cout << "## Initial Energy density: " << density << endl;
    
    random_gauge_trafo(U, gparams.seed);
    plaquette = gauge_energy(U);
    energy_density(U, density, Q);
    cout << "## Plaquette after rnd trafo: " << std::scientific << std::setw(15) << plaquette/U.getVolume()/double(U.getNc())/6. << endl; 
    cout << "## Energy density: " << density << endl;
    
    if(Wloop) {
      std::ostringstream os;
      os << "wilsonloop.";
      auto prevw = os.width(6);
      auto prevf = os.fill('0');
      os << i;
      os.width(prevw);
      os.fill(prevf);
      os << ".dat" << std::ends;
      compute_spacial_loops(U, os.str());
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
    if(lyapunov) {
      std::ostringstream os;
      os << "lyapunov.";
      auto prevw = os.width(6);
      auto prevf = os.fill('0');
      os << i;
      os.width(prevw);
      os.fill(prevf);

      // PRNG engine
      std::mt19937 engine(gparams.seed+i);
      gaugemonomial<double, su2> gm(0);
      kineticmonomial<double, su2> km(0);
      km.setmdpassive();
      
      std::list<monomial<double, su2>*> monomial_list;
      monomial_list.push_back(&gm);
      monomial_list.push_back(&km);
      md_params mdparams(n_steps, tau);
      integrator<double, su2> * md_integ = set_integrator<double, su2>(integs, exponent);

      compute_lyapunov(U, engine, mdparams, monomial_list, *md_integ, os.str(), exponent);
      delete(md_integ);
    }
  }

  return(0);
}
