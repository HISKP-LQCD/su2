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
#include<iomanip>
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

  for(size_t i = 500; i < 501; i+=20) {
    std::ostringstream os;
    os << "config." << Ls << "." << Lt << ".b" << U.getBeta() << "." << i << std::ends;
    U.load(os.str());
    
    double plaquette = gauge_energy(U);
    cout << "## Initital Plaquette: " << plaquette/U.getVolume()/N_c/6. << endl; 
    cout << "## Initial Energy density: " << energy_density(U) << endl;
    
    random_gauge_trafo(U, 654321);
    plaquette = gauge_energy(U);
    cout << "## Plaquette after rnd trafo: " << std::scientific << std::setw(15) << plaquette/U.getVolume()/N_c/6. << endl; 
    cout << "## Energy density: " << energy_density(U) << endl;
    
    {
      std::ostringstream os;
      os << "wilsonloop.";
      auto prevw = os.width(6);
      auto prevf = os.fill('0');
      os << i;
      os.width(prevw);
      os.fill(prevf);
      os << ".dat" << std::ends;
      compute_all_loops(U, os.str());
    }
    {
      std::ostringstream os;
      os << "gradient_flow.";
      auto prevw = os.width(6);
      auto prevf = os.fill('0');
      os << i;
      os.width(prevw);
      os.fill(prevf);
      gradient_flow(U, os.str(), 7.99);
    }
  }

  return(0);
}
