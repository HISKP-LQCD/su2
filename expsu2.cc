#include"expsu2.hh"
#include"su2.hh"
#include<vector>
#include<complex>

_su2 exp(std::vector<double> const & x) {
  const double alpha = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  std::vector<double> n = {x[0]/alpha, x[1]/alpha, x[2]/alpha};
  const double salpha = sin(alpha);
  _su2 res(Complex(cos(alpha), salpha*n[2]), 
           salpha*Complex(n[1], n[0]));
  return res;
}
