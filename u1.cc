#include"u1.hh"
#include<complex>


_u1 operator*(const _u1 &U1, const _u1 &U2) {
  _u1 res;
  res.a = U1.a + U2.a;
  return(res);
}

Complex operator*(const _u1 &U1, const Complex &U2) {
  return(std::exp(2*U1.geta()*M_PI*Complex(0., 1.)) * U2);
}
Complex operator*(const Complex &U1, const _u1 &U2) {
  return(U1 * std::exp(2*U2.geta()*M_PI*Complex(0., 1.)));
}


Complex operator+(const _u1 &U1, const _u1 &U2) {
  return(std::exp(2*U1.a*M_PI*Complex(0., 1.)) +
         std::exp(2*U2.a*M_PI*Complex(0., 1.)));
}

Complex operator-(const _u1 &U1, const _u1 &U2) {
  return(std::exp(2*U1.a*M_PI*Complex(0., 1.)) -
         std::exp(2*U2.a*M_PI*Complex(0., 1.)));
}

//Complex operator+=(const Complex &c, const _u1 &U) {
//  return c + std::exp(2*U.geta()*M_PI*Complex(0., 1.));
//}

void operator+=(Complex & U1, const _u1 & U2) {
  U1 += std::exp(2*U2.geta()*M_PI*Complex(0., 1.));
}

double trace(Complex c) {
  return(std::real(c));
}
