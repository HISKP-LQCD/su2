#pragma once 

#include<complex>
//#include<cmath>

using Complex = std::complex<double>;

const size_t N_c = 2;

class _su2 {
public:
  explicit _su2() : a(0), b(0) {}
  explicit _su2(Complex a, Complex b) : a(a), b(b) {}
  _su2(const _su2& U) : a(U.a), b(U.b) {}

  friend _su2 operator+(const _su2 &U1, const _su2 &U2);
  friend _su2 operator-(const _su2 &U1, const _su2 &U2);
  friend _su2 operator*(const _su2 &U1, const _su2 &U2);
  _su2& operator*=(const _su2 &U1) {
    Complex a = this->a;
    this->a = a*U1.a - this->b*std::conj(U1.b);
    this->b = a*U1.b + this->b*std::conj(U1.a);
    return *this;
  }
  _su2 round(size_t n) const {
    double dn = n;
    return _su2(Complex(std::round(std::real(a)*dn), std::round(std::imag(a)*dn))/dn,
                Complex(std::round(std::real(b)*dn), std::round(std::imag(b)*dn))/dn);
  }

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
  void operator+=(const _su2 &U) {
    a += U.a;
    b += U.b;
  }
  void set(const Complex _a, const Complex _b) {
    a = _a;
    b = _b;
  }
  _su2 dagger() const {
    return(_su2(std::conj(a), -b));
  }
  double trace() {
    return(2.*std::real(a));
  }
  Complex det() {
    return(a*std::conj(a) + b*std::conj(b)); 
  }
  void restoreSU() {
    double r = sqrt(std::abs(a)*std::abs(a) + std::abs(b)*std::abs(b));
    a /= r;
    b /= r;
  }

private:
  Complex a, b;
};

template<class matrix> double trace(matrix const &U) {
  double a = std::real(U.geta());
  return(2*a);
}

_su2 operator*(const _su2 &U1, const _su2 &U2);
_su2 operator+(const _su2 &U1, const _su2 &U2);
_su2 operator-(const _su2 &U1, const _su2 &U2);

using su2 = _su2;
