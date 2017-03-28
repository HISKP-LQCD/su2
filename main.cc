#include<iostream>
#include<vector>
#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"

using std::vector;
using std::cout;
using std::endl;


int main() {
  const size_t Ls = 8, Lt = 16;
  const double beta = 5.5;
  auto U = hotstart(Ls, Lt, 123456, 0.15);
  //auto config = coldstart(Ls, Lt);
  
  double plaquette = gauge_energy(U);
  cout << "Initital Plaquette: " << plaquette/U.getVolume()/2./6. << endl; 
  random_gauge_trafo(U, 654321);

  plaquette = gauge_energy(U);
  cout << "Plaquette after rnd trafo: " << plaquette/U.getVolume()/2./6. << endl; 
  return(0);
}

