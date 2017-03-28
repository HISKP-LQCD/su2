#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"

#include<iostream>
#include<vector>

using std::vector;
using std::cout;
using std::endl;


int main() {
  const size_t Ls = 8, Lt = 16;
  const double beta = 4.5;
  const size_t N_hit = 10;
  const size_t N_meas = 2000;
  const double delta = 0.1;
  auto U = hotstart(Ls, Lt, 123456, 0.2);
  //auto config = coldstart(Ls, Lt);
  
  double plaquette = gauge_energy(U);
  cout << "Initital Plaquette: " << plaquette/U.getVolume()/N_c/6. << endl; 

  random_gauge_trafo(U, 654321);
  plaquette = gauge_energy(U);
  cout << "Plaquette after rnd trafo: " << plaquette/U.getVolume()/N_c/6. << endl; 

  double rate = 0.;
  for(size_t i = 0; i < N_meas; i++) {
    rate += sweep(U, 13243546, delta, N_hit, beta);
    //cout << "Plaquette after sweep: " << i << " " << gauge_energy(U)/U.getVolume()/N_c/6. << endl;
    cout << i << " " << gauge_energy(U)/U.getVolume()/N_c/6. << endl;
  }
  cout << rate/static_cast<double>(N_meas) << endl;
  return(0);
}

