#pragma once

#include"su2.hh"
#include"u1.hh"
#include<vector>
#include<random>
#include<cassert>
#include<cmath>

template<typename Float> class adjointsu2 {
public:
  adjointsu2(Float _a, Float _b, Float _c) : a(_a), b(_b), c(_c) {}
  adjointsu2() : a(0.), b(0.), c(0.) {}
  void flipsign() {
    a = -a;
    b = -b;
    c = -c;
  }
  Float geta() const {
    return a;
  }
  Float getb() const {
    return b;
  }
  Float getc() const {
    return c;
  }
  void seta(Float _a) {
    a = _a;
  }
  void setb(Float _a) {
    b = _a;
  }
  void setc(Float _a) {
    c = _a;
  }
  void setzero() {
    a = b = c = 0.;
  }
  adjointsu2<Float> round(size_t n) const {
    Float dn = n;
    return adjointsu2(std::round(a * dn) / dn, std::round(b * dn) / dn, std::round(c * dn) / dn);
  }
  void operator=(const adjointsu2 &A) {
    a = A.geta();
    b = A.getb();
    c = A.getc();
  }
  void operator+=(const adjointsu2 &A) {
    a += A.geta();
    b += A.getb();
    c += A.getc();
  }
  void operator-=(const adjointsu2 &A) {
    a -= A.geta();
    b -= A.getb();
    c -= A.getc();
  }
  
private:
  Float a, b, c;
};

template<typename Float=double> inline adjointsu2<Float> get_deriv(su2 & A) {
  const Complex a = A.geta(), b = A.getb();
  return adjointsu2<Float>(2.*std::imag(b), 2.*std::real(b), 2.*std::imag(a));
}

template<typename Float> class adjointu1 {
public:
  adjointu1(Float _a) : a(_a) {}
  adjointu1() : a(0.) {}
  void flipsign() {
    a = -a;
  }
  Float geta() const {
    return a;
  }
  void seta(Float _a) {
    a = _a;
  }
  void setzero() {
    a = 0.;
  }

  adjointu1<Float> round(size_t n) const {
    Float dn = n;
    return adjointu1(std::round(a * dn) / dn);
  }
  void operator=(const adjointu1 &A) {
    a = A.geta();
  }
  void operator+=(const adjointu1 &A) {
    a += A.geta();
  }
  void operator-=(const adjointu1 &A) {
    a -= A.geta();
  }
  
private:
  Float a;
};

template<typename Float=double> inline adjointu1<Float> get_deriv(Complex & A) {
  return adjointu1<Float>(std::imag(A));
}

// The following class will be used to deliver the
// adjoint type depending on the gauge group
template<typename Float, class Group> struct adjoint_type {
  typedef Group type;
};

template<typename Float> struct adjoint_type<Float, su2> {
  typedef adjointsu2<Float> type;
};

template<typename Float> struct adjoint_type<Float, _u1> {
  typedef adjointu1<Float> type;
};


template<typename Float, class Group, typename lint=int> class adjointfield {
public:
  using value_type = typename adjoint_type<Float, Group>::type;
  adjointfield(const size_t Lx, const size_t Ly, const size_t Lz, const size_t Lt, const size_t ndims = 4) : 
    Lx(Lx), Ly(Ly), Lz(Lz), Lt(Lt), volume(Lx*Ly*Lz*Lt), ndims(ndims) {
    data.resize(volume*4);
  }
  adjointfield(const adjointfield &U) :
    Lx(U.getLx()), Ly(U.getLy()), Lz(U.getLz()), Lt(U.getLt()), volume(U.getVolume()), ndims(U.getndims()) {
    data.resize(volume*ndims);
    for(size_t i = 0; i < getSize(); i++) {
      data[i] = U[i];
    }
  }
  void flipsign() {
    for(size_t i = 0; i < getSize(); i++) {
      data[i].flipsign();
    }
  }
  size_t storage_size() const { return data.size() * sizeof(value_type); };
  size_t getLx() const {
    return(Lx);
  }
  size_t getLy() const {
    return(Ly);
  }
  size_t getLz() const {
    return(Lz);
  }
  size_t getLt() const {
    return(Lt);
  }
  size_t getndims() const {
    return(ndims);
  }
  size_t getVolume() const {
    return(volume);
  }
  size_t getSize() const {
    return(volume*ndims);
  }
  void operator=(const adjointfield<Float, Group> &U) {
    Lx = U.getLx();
    Ly = U.getLz();
    Lz = U.getLz();
    Lt = U.getLt();
    volume = U.getVolume();
    data.resize(U.getSize());
    for(size_t i = 0; i < U.getSize(); i++) {
      data[i] = U[i];
    }
  }

  value_type &operator()(lint const t, lint const x, lint const y, lint const z, size_t const mu) {
    return data[ getIndex(t, x, y, z, mu) ];
  }

  const value_type &operator()(lint const t, lint const x, lint const y, lint const z, size_t const mu) const {
    return data[ getIndex(t, x, y, z, mu) ];
  }

  value_type &operator()(std::vector<lint> const &coords, size_t const mu) {
    return data[ getIndex(coords[0], coords[1], coords[2], coords[3], mu) ];
  }

  const value_type &operator()(std::vector<lint> const &coords, size_t const mu) const {
    return data[ getIndex(coords[0], coords[1], coords[2], coords[3], mu) ];
  }

  value_type &operator[](size_t const index) {
    return data[ index ];
  }

  const value_type &operator[](size_t const index) const {
    return data[ index ];
  }

private:
  size_t Lx, Ly, Lz, Lt, volume, ndims;
  
  std::vector<value_type> data;

  size_t getIndex(const lint t, const lint x, const lint y, const lint z, const size_t mu) const {
    size_t y0 = (t + Lt) % Lt;
    size_t y1 = (x + Lx) % Lx;
    size_t y2 = (y + Ly) % Ly;
    size_t y3 = (z + Lz) % Lz;
    size_t _mu = (mu + ndims) % ndims;
    return( (((y0*Lx + y1)*Ly + y2)*Lz + y3)*ndims + _mu );
  }
};


template<typename Float> inline adjointsu2<Float> operator*(const Float &x, const adjointsu2<Float> &A) {
  return adjointsu2<Float>(x * A.geta(), x * A.getb(), x * A.getc());
}

template<typename Float> inline adjointu1<Float> operator*(const Float &x, const adjointu1<Float> &A) {
  return adjointu1<Float>(x * A.geta());
}


template<typename Float, class Group>  adjointfield<Float, su2> operator*(const Float &x, const adjointfield<Float, Group> &A) {
  adjointfield<Float, Group> res(A.getLx(), A.getLy(), A.getLz(), A.getLt());
  for(size_t i = 0; i < A.getSize(); i++) {
    res[i] = x * A[i];
  }
  return res;
}

template<typename Float> Float operator*(const adjointfield<Float, su2> &A, const adjointfield<Float, su2> &B) {
  Float res = 0.;
  assert(A.getSize() == B.getSize());
  for(size_t i = 0; i < A.getSize(); i++) {
    res += A[i].geta()*B[i].geta() + A[i].getb()*B[i].getb() + A[i].getc()*B[i].getc();
  }
  return res;
}

template<typename Float> Float operator*(const adjointfield<Float, _u1> &A, const adjointfield<Float, _u1> &B) {
  Float res = 0.;
  assert(A.getSize() == B.getSize());
  for(size_t i = 0; i < A.getSize(); i++) {
    res += A[i].geta()*B[i].geta();
  }
  return res;
}


template<class URNG, typename Float> void initnormal(URNG &engine, adjointfield<Float, su2> &A) {
  std::normal_distribution<double> normal(0., 1.);
  for(size_t i = 0; i < A.getSize(); i++) {
    A[i].seta(Float(normal(engine)));
    A[i].setb(Float(normal(engine)));
    A[i].setc(Float(normal(engine)));
  }
  return;
}

template<class URNG, typename Float> void initnormal(URNG &engine, adjointfield<Float, _u1> &A) {
  std::normal_distribution<double> normal(0., 1.);
  for(size_t i = 0; i < A.getSize(); i++) {
    A[i].seta(Float(normal(engine)));
  }
  return;
}

template<typename Float, class Group> inline void zeroadjointfield(adjointfield<Float, Group> &A) {
  for(size_t i = 0; i < A.getSize(); i++) {
    A[i].setzero();
  }
}

