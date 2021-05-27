#pragma once 

#include"accum_type.hh"
#include"traceless_antiherm.hh"
#include<complex>


template<typename Float=double> class _su2 {
public:
  typedef std::complex<Float> Complex;
  const size_t N_c = 2;
  explicit _su2() : a(0), b(0) {}
  explicit _su2(Complex a, Complex b) : a(a), b(b) {}
  _su2(const _su2<Float>& U) : a(U.a), b(U.b) {}

  template<typename F>
  friend inline _su2<F> operator+(const _su2<F> &U1, const _su2<F> &U2);
  template<typename F>
  friend inline _su2<F> operator-(const _su2<F> &U1, const _su2<F> &U2);
  template<typename F>
  friend inline _su2<F> operator*(const _su2<F> &U1, const _su2<F> &U2);
  _su2<Float>& operator*=(const _su2<Float> &U1) {
    Complex a = this->a;
    this->a = a*U1.a - this->b*std::conj(U1.b);
    this->b = a*U1.b + this->b*std::conj(U1.a);
    return *this;
  }
  _su2<Float> round(size_t n) const {
    double dn = n;
    return _su2<Float>(Complex(std::round(std::real(a)*dn), std::round(std::imag(a)*dn))/dn,
                       Complex(std::round(std::real(b)*dn), std::round(std::imag(b)*dn))/dn);
  }

  inline Complex geta() const {
    return(a);
  }
  inline Complex getb() const {
    return(b);
  }
  inline void operator=(const _su2<Float> &U) {
    a = U.geta();
    b = U.getb();
  }
  inline void operator+=(const _su2<Float> &U) {
    a += U.a;
    b += U.b;
  }
  void set(const Complex _a, const Complex _b) {
    a = _a;
    b = _b;
  }
  inline _su2<Float> dagger() const {
    return(_su2<Float>(std::conj(a), -b));
  }
  inline double retrace() {
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

template<typename Float=double>
inline double retrace(_su2<Float> const &U) {
  double a = std::real(U.geta());
  return(2*a);
}
template<typename Float=double>
inline typename _su2<Float>::Complex trace(_su2<Float> const &U) {
  double a = std::real(U.geta());
  return(typename _su2<Float>::Complex(2*a, 0.));
}

template<typename Float=double>
inline _su2<Float> traceless_antiherm(const _su2<Float>& x) {
  return(_su2<Float>(0.5*(x.geta()-std::conj(x.geta())), x.getb()));
}

template<typename Float=double>
inline _su2<Float> operator*(const _su2<Float> &U1, const _su2<Float> &U2) {
  _su2<Float> res;
  res.a = U1.a*U2.a - U1.b*std::conj(U2.b);
  res.b = U1.a*U2.b + U1.b*std::conj(U2.a);
  return(res);
}

template<typename Float=double>
inline _su2<Float> operator+(const _su2<Float> &U1, const _su2<Float> &U2) {
  _su2<Float> res;
  res.a = U1.a + U2.a;
  res.b = U1.b + U2.b;
  return(res);
}

template<typename Float=double>
inline _su2<Float> operator-(const _su2<Float> &U1, const _su2<Float> &U2) {
  _su2<Float> res;
  res.a = U1.a - U2.a;
  res.b = U1.b - U2.b;
  return(res);
}


using su2 = _su2<double>;

