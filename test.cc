#include"su2.hh"
#include"geometry.hh"
#include<iostream>
#include<vector>

using std::vector;
using std::cout;
using std::endl;


int main() {
  size_t Ls = 8, Lt = 16;

  
  vector<su2> config;
  config.resize(Ls*Ls*Ls*Lt);
  
  // set all to 1
  for(auto i = config.begin(), end =config.end(); i < end; i++) {
    *i = su2(1., 0.);
  }
  
  cout << config[0].trace() << " " << config[0].det() << " " << config[0].geta() << " " << config[0].getb() << endl;

  su2 U1;
  U1 = config[0] * config[1];
  cout << U1.det() << endl;

  U1 = su2(Complex(0.8, 0.3), Complex(0.1, 0.4));
  su2 U2(0.9, Complex(0.3, 0.2));
  U1.rescale();
  U2.rescale();
  su2  U3 = U1 * U2;
  cout << U3.det() << " " << trace<su2>(U1 * U2) << " " << U3.geta() << " " << U3.getb() << endl;

  geometry geo(Ls, Lt);
  size_t i = geo.getIndex(7, 5, 3, 9);
  size_t c[4];
  geo.getCoordinate(c, i);
  cout << i << " " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << endl;

  su2 U = config[0] * config[0].dagger();
  cout << U.geta() << " " << U.getb() << endl;

  return(0);
}
