#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"wilsonloop.hh"
#include"md_update.hh"
#include"monomial.hh"
#include"gradient_flow.hh"

#include<iostream>
#include<vector>
#include<random>

using std::vector;
using std::cout;
using std::endl;


int main() {
  const size_t Ls = 8, Lt = 16;
  const double beta = 4.5;
  const size_t N_meas = 10;
  const size_t N_save = 1;
  const size_t N_rev = 100;
  const int seed = 13526463;
  gaugeconfig U(Ls, Lt, beta);
  U = hotstart(Ls, Lt, 123456, 0.10);

  md_params params(100, 1.0);
  
  std::mt19937 engine(seed);

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

  double rate = 0.;
  for(size_t i = 0; i < N_meas; i++) {
    params.disablerevtest();
    if((i+1) % N_rev == 0) {
      params.enablerevtest();
    }
    md_update(U, engine, params, monomial_list);

    rate += params.getaccept();
    cout << i+1 << " " << gauge_energy(U)/U.getVolume()/N_c/6. << " " << params.getdeltaH() << " " 
         << params.getaccept() << " " << rate/static_cast<double>(i+1) << " ";
    if(params.getrevtest()) {
      cout << params.getdeltadeltaH();
    }
    else cout << "NA";
    cout << endl;

    if(i > 0 && i % N_save == 0) {
      {
        std::ostringstream os;
        os << "wilsonloop." << i << ".dat" << std::ends;
        compute_all_loops(U, os.str());
      }
      {
        std::ostringstream os;
        os << "gradient_flow." << i << ".dat" << std::ends;
        gradient_flow(U, os.str());
      }
    }
  }
  cout << rate/static_cast<double>(N_meas) << endl;
  return(0);
}
