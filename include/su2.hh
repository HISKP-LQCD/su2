#pragma once 

#include<complex>
#include<cmath>

using Complex = std::complex<double>;

class _su2 {
public:
  explicit _su2() : a(0), b(0) {}
  explicit _su2(Complex a, Complex b) : a(a), b(b) {}
  _su2(const _su2& U) : a(U.a), b(U.b) {}

  friend _su2 operator+(const _su2 &U1, const _su2 &U2);
  friend _su2 operator-(const _su2 &U1, const _su2 &U2);
  friend _su2 operator*(const _su2 &U1, const _su2 &U2);

  Complex geta() const {
    return(a);
  }
  Complex getb() const {
    return(b);
  }
  void operator=(const _su2 &U) {
    a = U.geta();
    b = U.getb();
  }
  void set(const Complex _a, const Complex _b) {
    a = _a;
    b = _b;
  }
  _su2 dagger() const {
    return(_su2(a, b));
  }
  Complex trace() {
    return(2.*std::real(a));
  }
  Complex det() {
    return(a*std::conj(a) + b*std::conj(b)); 
  }
  void rescale() {
    double r = sqrt(std::abs(a)*std::abs(a) + std::abs(b)*std::abs(b));
    a /= r;
    b /= r;
  }

private:
  Complex a, b;
};

Complex trace(const _su2 &U) {
  double a = std::real(U.geta());
  return(2*a);
}

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

using su2 = _su2;
