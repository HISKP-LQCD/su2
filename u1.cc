#include"u1.hh"
#include<complex>


_u1 operator*(const _u1 &U1, const _u1 &U2) {
  _u1 res;
  res.a = U1.a + U2.a;
  return(res);
}

Complex operator+(const _u1 &U1, const _u1 &U2) {
  return(std::exp(2*U1.a*M_PI*Complex(0., 1.)) +
         std::exp(2*U2.a*M_PI*Complex(0., 1.)));
}

Complex operator-(const _u1 &U1, const _u1 &U2) {
  return(std::exp(2*U1.a*M_PI*Complex(0., 1.)) -
         std::exp(2*U2.a*M_PI*Complex(0., 1.)));
}

