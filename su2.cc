#include"su2.hh"
#include<complex>


_su2 operator*(const _su2 &U1, const _su2 &U2) {
  _su2 res;
  res.a = U1.a*U2.a - U1.b*std::conj(U2.b);
  res.b = U1.a*U2.b + U1.b*std::conj(U2.a);
  return(res);
}

_su2 operator+(const _su2 &U1, const _su2 &U2) {
  _su2 res;
  res.a = U1.a + U2.a;
  res.b = U1.b + U2.b;
  return(res);
}

_su2 operator-(const _su2 &U1, const _su2 &U2) {
  _su2 res;
  res.a = U1.a - U2.a;
  res.b = U1.b - U2.b;
  return(res);
}

