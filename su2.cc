#include"su2.hh"
#include<complex>

double trace(const su2 &U) {
  double a = std::real(U.geta());
  return(2*a);
}

su2 operator*(const su2 &U1, const su2 &U2) {
  su2 res;
  res.a = U1.a*U2.a - U1.b*std::conj(U2.b);
  res.b = U1.a*U2.b + U1.b*std::conj(U2.a);
  return(res);
}

su2 operator+(const su2 &U1, const su2 &U2) {
  su2 res;
  res.a = U1.a + U2.a;
  res.b = U1.b + U2.b;
  return(res);
}

su2 operator-(const su2 &U1, const su2 &U2) {
  su2 res;
  res.a = U1.a - U2.a;
  res.b = U1.b - U2.b;
  return(res);
}
