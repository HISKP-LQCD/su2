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
  const double beta = 4.5;
  const size_t N_meas = 500;
  const size_t N_save = 20;
  const size_t N_rev = 1;
  const int seed = 13526463;
  gaugeconfig U(Ls, Lt, beta);
  U = hotstart(Ls, Lt, 123456, 0.10);

  md_params params(1, 0.02);
  
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
