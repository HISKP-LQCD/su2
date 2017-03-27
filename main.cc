#include<iostream>
#include<vector>
#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"

using std::vector;
using std::cout;
using std::endl;


int main() {
  const size_t Ls = 8, Lt = 16;
  const double beta = 5.5;
  auto config = hotstart(Ls, Lt, 123456, 0.5);
  //auto config = coldstart(Ls, Lt);
  
  double plaquette = gauge_energy(config);

  cout << "Initital Plaquette: " << plaquette/config.getVolume()/2./6. << endl; 
  return(0);
}

