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

#include<iostream>
#include<sstream>
#include<vector>
#include<random>

using std::vector;
using std::cout;
using std::endl;


int main() {
  const size_t Ls = 8, Lt = 8;
  const double beta = 2.3;
  const size_t N_meas = 500;
  const size_t N_save = 20;
  const size_t N_rev = 1;
  const int seed = 13526463;
  gaugeconfig U(Ls, Lt, beta);
  U = hotstart(Ls, Lt, 123456, 0.10);

  md_params params(50, 1.0);
  
  std::mt19937 engine(seed);

  std::ostringstream os;
  os << "config." << Ls << "." << Lt << "." << U.getBeta() << ".1980" << std::ends;
  //os << "gauge-links-0000.bin" << std::ends;
  //os << "indexing.bin" << std::ends;
  U.load(os.str());

  double plaquette = gauge_energy(U);
  cout << "## Initital Plaquette: " << plaquette/U.getVolume()/N_c/6. << endl; 
  cout << "## Initial Energy density: " << energy_density(U) << endl;

  random_gauge_trafo(U, 654321);
  plaquette = gauge_energy(U);
  cout << "## Plaquette after rnd trafo: " << plaquette/U.getVolume()/N_c/6. << endl; 
  cout << "## Energy density: " << energy_density(U) << endl;

  size_t i = 0;
  {
    std::ostringstream os;
    os << "wilsonloop." << i+1 << ".dat" << std::ends;
    compute_all_loops(U, os.str());
  }
  {
    std::ostringstream os;
    os << "gradient_flow." << i+1 << ".dat" << std::ends;
    gradient_flow(U, os.str());
  }

  return(0);
}
