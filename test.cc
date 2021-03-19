#include"su2.hh"
#include"u1.hh"
#include"random_element.hh"
#include<iostream>
#include<vector>

using std::vector;
using std::cout;
using std::endl;


int main() {
  size_t Ls = 8, Lt = 16;


  cout << "Tests of SU(2)" << endl << endl;;
  
  vector<su2> config;
  config.resize(Ls*Ls*Ls*Lt);
  
  // set all to 1
  for(auto i = config.begin(), end =config.end(); i < end; i++) {
    *i = su2(1., 0.);
  }

  su2 U(1., 0.);
  cout << "retrace, det, a and b of unit matrix" << endl;
  cout << U.retrace() << " = " << retrace(U) << endl << U.det()
       << " " << U.geta() << " " << U.getb() << endl;

  su2 U1;
  U1 = config[0] * config[1];
  cout << "det of U * U =" << U1.det() << endl;

  U1 = su2(Complex(0.8, 0.3), Complex(0.1, 0.4));
  su2 U2(0.9, Complex(0.3, 0.2));
  U1.restoreSU();
  U2.restoreSU();
  su2  U3 = U1 * U2;
  cout << "U3 = U1 * U2: det, retrace, a, b" << endl;
  cout << "should be: 1, 1.3264 (0.6632,0.184826) (0.293548,0.6632)" << endl;
  cout << U3.det() << " " << retrace(U1 * U2) << " " << U3.geta() << " " << U3.getb() << endl;

  U = config[0] * config[0].dagger();
  cout << "Test of U*U^dagger, should be: (1,0) (0,0)" << endl;
  cout << U.geta() << " " << U.getb() << endl;

  cout << endl << "Tests of U(1)" << endl << endl;;
  _u1 u;
  double a = u.retrace();
  cout << "test basic retrace, det and geta functions" << endl;
  cout << u.retrace() << " = " << retrace(u) << endl
       << u.det() << " " << u.geta() << " " << endl;

  _u1 x(0.5), y(0.7), z;
  cout << "test of det, the two following must be equal" << endl;
  cout << x.det() << " = " << std::exp(0.5*2*M_PI*Complex(0., 1.)) << endl;
  cout << "test of retrace, the two following must be equal" << endl;
  cout << x.retrace() << " = " << retrace(x) << endl;
  cout << "test multiplication, the two complex numbers must agree" << endl;

  z = x * y;
  cout << z.det() << " = " << x.det() * y.det() << endl;
  return(0);
}
