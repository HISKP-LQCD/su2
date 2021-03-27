#include"expsu2.hh"
#include"su2.hh"
#include"u1.hh"
#include"adjointfield.hh"
#include<vector>
#include<complex>

_su2 exp(adjointsu2<double> const & x) {
  double a = x.geta(), b = x.getb(), c = x.getc();
  // normalise the vector
  const double alpha = sqrt(a*a+b*b+c*c);
  std::vector<double> n = {a/alpha, b/alpha, c/alpha};

  const double salpha = sin(alpha);
  _su2 res(Complex(cos(alpha), salpha*n[2]), 
           salpha*Complex(n[1], n[0]));
  return res;
}

inline _u1 exp(adjointu1<double> const & x) {
  return _u1(x.geta());
}
