#pragma once

#include<vector>
#include<random>
#include<cassert>

template<class T> class adjoint {
public:
  adjoint(T _a, T _b, T _c) : a(_a), b(_b), c(_c) {}
  adjoint() : a(0.), b(0.), c(0.) {}
  void flipsign() {
    a = -a;
    b = -b;
    c = -c;
  }
  T geta() const {
    return a;
  }
  T getb() const {
    return b;
  }
  T getc() const {
    return c;
  }
  void seta(T _a) {
    a = _a;
  }
  void setb(T _a) {
    b = _a;
  }
  void setc(T _a) {
    c = _a;
  }
  void operator=(const adjoint &A) {
    a = A.geta();
    b = A.getb();
    c = A.getc();
  }
  void operator+=(const adjoint &A) {
    a += A.geta();
    b += A.getb();
    c += A.getc();
  }
  void operator-=(const adjoint &A) {
    a -= A.geta();
    b -= A.getb();
    c -= A.getc();
  }
  
private:
  T a, b, c;
};

template<class T> class adjointfield {
public:
  using value_type = adjoint<T>;
  
  adjointfield(const size_t Ls, const size_t Lt) : 
    Ls(Ls), Lt(Lt), volume(Ls*Ls*Ls*Lt) {
    data.resize(volume*4);
  }
  adjointfield(const adjointfield &U) :
    Ls(U.getLs()), Lt(U.getLt()), volume(U.getVolume()) {
    data.resize(volume*4);
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
  size_t getLs() const {
    return(Ls);
  }
  size_t getLt() const {
    return(Lt);
  }
  size_t getVolume() const {
    return(volume);
  }
  size_t getSize() const {
    return(volume*4);
  }
  void operator=(const adjointfield &U) {
    Ls = U.getLs();
    Lt = U.getLt();
    volume = U.getVolume();
    data.resize(U.getSize());
    for(size_t i = 0; i < U.getSize(); i++) {
      data[i] = U[i];
    }
  }

  value_type &operator()(size_t const t, size_t const x, size_t const y, size_t const z, size_t const mu) {
    return data[ getIndex(t, x, y, z, mu) ];
  }

  const value_type &operator()(size_t const t, size_t const x, size_t const y, size_t const z, size_t const mu) const {
    return data[ getIndex(t, x, y, z, mu) ];
  }

  value_type &operator()(std::vector<size_t> const &coords, size_t const mu) {
    return data[ getIndex(coords[0], coords[1], coords[2], coords[3], mu) ];
  }

  const value_type &operator()(std::vector<size_t> const &coords, size_t const mu) const {
    return data[ getIndex(coords[0], coords[1], coords[2], coords[3], mu) ];
  }

  value_type &operator[](size_t const index) {
    return data[ index ];
  }

  const value_type &operator[](size_t const index) const {
    return data[ index ];
  }

private:
  size_t Ls, Lt, volume;
  
  std::vector<value_type> data;
  
  size_t getIndex(const size_t t, const size_t x, const size_t y, const size_t z, const size_t mu) const {
    size_t y0 = (t + Lt) % Lt;
    size_t y1 = (x + Ls) % Ls;
    size_t y2 = (y + Ls) % Ls;
    size_t y3 = (z + Ls) % Ls;
    size_t _mu = (mu + 4) % 4;
    return( (((y0*Ls + y1)*Ls + y2)*Ls + y3)*4 + _mu );
  }
};


template<class T>  adjoint<T> operator*(const T &x, const adjoint<T> &A) {
  adjoint<T> res;
  res.seta(x * A.geta());
  res.setb(x * A.getb());
  res.setc(x * A.getc());
  return res;
}

template<class T>  adjointfield<T> operator*(const T &x, const adjointfield<T> &A) {
  adjointfield<T> res(A.getLs(), A.getLt());
  for(size_t i = 0; i < A.getSize(); i++) {
    res[i].seta( x * A[i].geta());
    res[i].setb( x * A[i].getb());
    res[i].setc( x * A[i].getc());
  }
  return res;
}

template<class T> T operator*(const adjointfield<T> &A, const adjointfield<T> &B) {
  T res = 0.;
  assert(A.getSize() == B.getSize());
  for(size_t i = 0; i < A.getSize(); i++) {
    res += A[i].geta()*B[i].geta() + A[i].getb()*B[i].getb() + A[i].getc()*B[i].getc();
  }
  return res;
}

template<class URNG, class T> adjointfield<T> initnormal(URNG &engine, size_t Ls, size_t Lt) {
  adjointfield<T> A(Ls, Lt);
  std::normal_distribution<double> normal(0., 1.);
  for(size_t i = 0; i < A.getSize(); i++) {
    A[i].seta(T(normal(engine)));
    A[i].setb(T(normal(engine)));
    A[i].setc(T(normal(engine)));
  }
  return A;
}

template<class T> void zeroadjointfield(adjointfield<T> &A) {
  for(size_t i = 0; i < A.getSize(); i++) {
    A[i].seta(0.);
    A[i].setb(0.);
    A[i].setc(0.);
  }
}
