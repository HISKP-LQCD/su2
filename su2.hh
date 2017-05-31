#pragma once 

#include<complex>

using Complex = std::complex<double>;

const size_t N_c = 2;

template<class E> class su2expression;

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

  Complex geta() const {
    return(a);
  }
  Complex getb() const {
    return(b);
  }
  void seta(Complex const &_a) {
    a = _a;
  }
  void setb(Complex const &_b) {
    b = _b;
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
  void rescale() {
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

// marker template class
template<class E> class su2expression {
public:
  E expr;

  explicit su2expression() : expr() {}
  explicit su2expression(E e) : expr(e) {}
  explicit su2expression(Complex a, Complex b): expr(a, b) {}
  Complex geta() const {
    return(expr.geta());
  }
  Complex getb() const {
    return(expr.getb());
  }
  void seta(Complex const &a) {
    expr.seta(a);
  }
  void setb(Complex const &b) {
    expr.setb(b);
  }
  template<class T> void operator=(const su2expression<T> &U) {
    expr.seta( U.geta() );
    expr.setb( U.getb() );
  }
  template<class T> void operator+=(const su2expression<T> &U) {
    expr.seta( expr.geta() += U.geta() );
    expr.setb( expr.getb() += U.getb() );
  }
  void set(const Complex _a, const Complex _b) {
    expr.seta(_a);
    expr.setb(_b);
  }
  su2expression<E> dagger() const {
    su2expression<E> tmp(std::conj(expr.geta()), -expr.getb());
    return(tmp);
  }
  double trace() {
    return(2.*std::real(expr.geta()));
  }
  Complex det() {
    return(expr.geta()*std::conj(expr.geta()) + expr.getb()*std::conj(expr.getb())); 
  }
  void rescale() {
    double r = sqrt(std::abs(expr.geta())*std::abs(expr.geta()) + std::abs(expr.getb())*std::abs(expr.getb()));
    expr.seta(expr.geta() /= r);
    expr.setb(expr.getb() /= r);
  }
};

// template class for add operations
template<class LE, class RE> class Add {
public:
  explicit Add(LE le, RE re) : lexpr(le), rexpr(re) {}
  Complex geta() const {
    return( lexpr.geta() + rexpr.geta() );
  }
  Complex getb() const {
    return( lexpr.getb() + rexpr.getb() );
  }
private:
  LE lexpr;
  RE rexpr;
};
// overload operator+
template<class LE, class RE> su2expression<Add<LE, RE>> operator+(su2expression<LE> le, su2expression<RE> re) {
  return su2expression<Add<LE, RE>>(Add<LE, RE>(le.expr, re.expr));
}

// template class for multiplications
template<class LE, class RE> class Mul {
public:
  explicit Mul(LE le, RE re) : lexpr(le), rexpr(re) {}
  Complex geta() const {
    return( lexpr.geta()*rexpr.geta()-lexpr.getb()*std::conj(rexpr.getb()) );
  }
  Complex getb() const {
    return( lexpr.geta()*rexpr.getb() + lexpr.getb()*std::conj(rexpr.geta()) );
  }
private:
  LE lexpr;
  RE rexpr;
};

// overload operator*
template<class LE, class RE> su2expression<Mul<LE, RE>> operator*(su2expression<LE> le, su2expression<RE> re) {
  return su2expression<Mul<LE, RE>>(Mul<LE, RE>(le.expr, re.expr));
}


using su2 = su2expression<_su2>;
