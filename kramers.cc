#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"wilsonloop.hh"
#include"kramers_md_update.hh"
#include"monomial.hh"
#include"gradient_flow.hh"
#include"energy_density.hh"

#include<iostream>
#include<sstream>
#include<vector>
#include<random>
#include <boost/program_options.hpp>

using std::vector;
using std::cout;
using std::endl;
using std::cerr;
namespace po = boost::program_options;

int main(int ac, char* av[]) {

  size_t Ls = 8, Lt = 16;
  double beta = 2.3;
  size_t N_meas = 2000;
  size_t N_rev = 1;
  double heat = 0.;
  size_t N_save = 200;
  int seed = 13526463;
  const size_t n_steps = 1;
  size_t exponent = 16;
  double tau = 1.;

  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce this help message")
      ("spatialsize,L", po::value<size_t>(&Ls), "spatial lattice size")
      ("temporalsize,T", po::value<size_t>(&Lt), "temporal lattice size")
      ("nsave", po::value<size_t>(&N_save)->default_value(1000), "N_save")
      //      ("nrev", po::value<size_t>(&N_rev)->default_value(1000), "N_rev")
      ("tau", po::value<double>(&tau)->default_value(1.), "trajectory length tau")
      ("nmeas,n", po::value<size_t>(&N_meas)->default_value(10), "total number of measurements")
      ("exponent", po::value<size_t>(&exponent)->default_value(16), "exponent for rounding")
      ("beta,b", po::value<double>(&beta), "beta value")
      ("seed,s", po::value<int>(&seed)->default_value(13526463), "PRNG seed")
      ("heat", po::value<double>(&heat)->default_value(1.), "randomness of the initial config, 1: hot, 0: cold")
      ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
      cout << desc << endl;
      return 1;
    }
    if (!vm.count("spatialsize") && !vm.count("help")) {
      std::cerr << "spatial lattice size must be given!" << endl;
      cout << endl << desc << endl;
      return 1;
    }
    if (!vm.count("temporalsize") && !vm.count("help")) {
      std::cerr << "temporal lattice size must be given!" << endl;
      cout << endl << desc << endl;
      return 1;
    }
    if (!vm.count("beta") && !vm.count("help")) {
      std::cerr << "beta value must be specified!" << endl;
      cout << endl << desc << endl;
      return 1;
    }
  }
  catch(std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  catch(...) {
    std::cerr << "Exception of unknown type!\n";
  }


  gaugeconfig U(Ls, Lt, beta);
  U = hotstart(Ls, Lt, seed, heat);

  md_params params(n_steps, tau);
  
  std::mt19937 engine(seed);

  double plaquette = gauge_energy(U);
  cout << "## Initital Plaquette: " << plaquette/U.getVolume()/N_c/6. << endl; 
  cout << "## Initial Energy density: " << energy_density(U) << endl;

  random_gauge_trafo(U, 654321);
  plaquette = gauge_energy(U);
  cout << "## Plaquette after rnd trafo: " << plaquette/U.getVolume()/N_c/6. << endl; 
  cout << "## Energy density: " << energy_density(U) << endl;

  // generate list of monomials
  gaugemonomial<double> gm(0);
  kineticmonomial<double> km(0);
  km.setmdpassive();

  std::list<monomial<double>*> monomial_list;
  monomial_list.push_back(&gm);
  monomial_list.push_back(&km);

  params.setkmax(5);

  double rate = 0.;
  for(size_t i = 0; i < N_meas; i++) {
    params.disablerevtest();
    kramers_md_update(U, engine, params, monomial_list);

    rate += params.getaccept();
    std::cout << i << " " << params.getaccept() << " " << gauge_energy(U)/U.getVolume()/N_c/6. << " " << params.getdeltaH() << " " 
              << rate/static_cast<double>(i+1) << std::endl;

    if(i > 0 && (i % N_save) == 0) {
      std::ostringstream os;
      os << "config." << Ls << "." << Lt << ".b" << beta << "." << i << std::ends;
      U.save(os.str());
    }
  }
  cout << "## Acceptance rate: " << rate/static_cast<double>(N_meas) << endl;

  std::ostringstream os;
  os << "config." << Ls << "." << Lt << ".b" << U.getBeta() << ".final" << std::ends;
  U.save(os.str());
  return(0);
}
