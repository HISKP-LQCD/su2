#pragma once 

#include"accum_type.hh"
#include"traceless_antiherm.hh"
#include<complex>
#include<iostream>
#include<fstream>
//#include<cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


using Complex = std::complex<double>;

//const size_t N_c = 2;

class _u1 {
public:
  const size_t N_c = 1;
  explicit _u1() : a(0) {}
  explicit _u1(double a) : a(a) {}
  _u1(const _u1& U) : a(U.a) {}

  friend Complex operator+(const _u1 &U1, const _u1 &U2);
  friend Complex operator-(const _u1 &U1, const _u1 &U2);
  friend _u1 operator*(const _u1 &U1, const _u1 &U2);
  // implicit conversion operator to complex
  operator Complex() const {
    return(std::exp(2*a*M_PI*Complex(0., 1.))); 
  }
  _u1& operator*=(const _u1 &U1) {
    this->a += U1.a;
    return *this;
  }

  double geta() const {
    return(a);
  }
  void operator=(const _u1 &U) {
    a = U.a;
  }
  void set(const double _a) {
    a = _a;
  }
  _u1 dagger() const {
    return(_u1(-a));
  }
  double retrace() {
    return(std::cos(a*2*M_PI));
  }
  Complex det() {
    return(std::exp(2*a*M_PI*Complex(0., 1.))); 
  }
  void restoreSU() {
  }

private:
  double a;
};

inline double retrace(_u1 const &U) {
  return(std::cos(U.geta()*2*M_PI));
}

inline double retrace(const Complex c) {
  return(std::real(c));
}

inline Complex trace(const Complex c) {
  return(c);
}

template<> struct accum_type<_u1> {
  typedef Complex type;
};

template<> inline Complex traceless_antiherm(const Complex& x) {
  return(Complex(0., std::imag(x)));
}

_u1 operator*(const _u1 &U1, const _u1 &U2);
Complex operator*(const _u1 &U1, const Complex &U2);
Complex operator*(const Complex &U1, const _u1 &U2);
Complex operator+(const _u1 &U1, const _u1 &U2);
Complex operator-(const _u1 &U1, const _u1 &U2);
void operator+=(Complex & U1, const _u1 & U2);

