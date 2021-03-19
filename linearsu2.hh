#pragma once

#include<complex>
#include<math.h>
#include<cmath>
#include<cassert>
#include"su2.hh"

using Complex = std::complex<double>;

class _Lsu2 {
public:
  const size_t N_c = 2;
  explicit _Lsu2() : m(0) {
    j[0] = m;
    j[1] = 0;
    j[2] = 0;
    j[3] = 0;
    s[0] = +1;
    s[1] = +1;
    s[2] = +1;
    s[3] = +1;
  }
  explicit _Lsu2(size_t m) : m(m) {
    j[0] = m;
    j[1] = 0;
    j[2] = 0;
    j[3] = 0;
    s[0] = +1;
    s[1] = +1;
    s[2] = +1;
    s[3] = +1;
    assert((m == j[0] + j[1] + j[2] + j[3]));
    assert(((abs(s[0]) == 1) && (abs(s[1]) == 1) && (abs(s[2]) == 1) && (abs(s[3]) == 1)));
  }
  explicit _Lsu2(size_t _m,
                 size_t * _j,
                 int * _s){
    for(int i = 0; i < 4; i++) {
      j[i] = _j[i];
      s[i] = _s[i];
    }
    m=_m;
    assert((m == (j[0] + j[1] + j[2] + j[3])));
    assert(((abs(s[0]) == 1) && (abs(s[1]) == 1) && (abs(s[2]) == 1) && (abs(s[3]) == 1)));
  }
  explicit _Lsu2(size_t m,
                 size_t j1, size_t j2, size_t j3,
                 int s1=1, int s2=1, int s3=1, int s4=1) : m(m) {
    j[0] = j1;
    j[1] = j2;
    j[2] = j3;
    s[0] = s1;
    s[1] = s2;
    s[2] = s3;
    s[3] = s4;
    j[3] = m - j[0] - j[1] - j[2];
    assert((m == (j[0] + j[1] + j[2] + j[3])));
    assert(((abs(s[0]) == 1) && (abs(s[1]) == 1) && (abs(s[2]) == 1) && (abs(s[3]) == 1)));
  }
  friend su2 operator*(const _Lsu2 &U1, const _Lsu2 &U2);
  friend su2 operator*(const su2 &U1, const _Lsu2 &U2);
  friend su2 operator*(const _Lsu2 &U1, const su2 &U2);

  size_t getj(const size_t i) const {
    return(j[i]);
  }
  void setjpair(const size_t i1, const size_t i2,
               const size_t _j1, const size_t _j2) {
    j[i1] = _j1;
    j[i2] = _j2;
    assert((m == (j[0] + j[1] + j[2] + j[3])));
    return;
  }
  int gets(const size_t i) const {
    return(s[i]);
  }
  void sets(const size_t i, const int _s) {
    s[i] = _s;
    return;
  }
  size_t getm() const {
    return(m);
  }
  void setm(const size_t _m) {
    m = _m;
  }
  double getJ() const {
    return(sqrt(static_cast<double>(j[0]*j[0] + j[1]*j[1] + j[2]*j[2] + j[3]*j[3])));
  }
  void operator=(const _Lsu2 &U) {
    for(int i = 0; i < 4; i++) {
      j[i] = U.getj(i);
      s[i] = U.gets(i);
      m = U.getm();
    }
    assert((m == (j[0] + j[1] + j[2] + j[3])));
    assert(((abs(s[0]) == 1) && (abs(s[1]) == 1) && (abs(s[2]) == 1) && (abs(s[3]) == 1)));
  }
  void set(const size_t * const _j, const int * const _s, const size_t _m) {
    for(int i = 0; i < 4; i++) {
      j[i] = _j[i];
      s[i] = _s[i];
    }
    m=_m;
    assert((m == (j[0] + j[1] + j[2] + j[3])));
    assert(((abs(s[0]) == 1) && (abs(s[1]) == 1) && (abs(s[2]) == 1) && (abs(s[3]) == 1)));
  }
  _Lsu2 dagger() const {
    return(_Lsu2(m, j[0], j[1], j[2], s[0], -s[1], -s[2], -s[3]));
  }
  double retrace() {
    return(2 * s[0] * static_cast<double>(j[0])/getJ());
  }
  Complex det() {
    return((static_cast<double>(j[0])*static_cast<double>(j[0]) +
            static_cast<double>(j[1])*static_cast<double>(j[1]) +
            static_cast<double>(j[2])*static_cast<double>(j[2]) +
            static_cast<double>(j[3])*static_cast<double>(j[3])) /
           static_cast<double>(getJ()*getJ())); 
  }
  void restoreSU() {
    assert((m == (j[0] + j[1] + j[2] + j[3])));
    assert(((abs(s[0]) == 1) && (abs(s[1]) == 1) && (abs(s[2]) == 1) && (abs(s[3]) == 1)));
    return;
  }
  double weight() {
    const double w = sqrt(2.)/getJ();
    return(w*w*w);
  }

private:
  size_t m;
  size_t j[4];
  int s[4];
};

template<> struct accum_type<_Lsu2> {
  typedef _su2 type;
};


su2 operator*(const _Lsu2 &U1, const _Lsu2 &U2);
su2 operator*(const su2 &U1, const _Lsu2 &U2);
su2 operator*(const _Lsu2 &U1, const su2 &U2);

using Lsu2 = _Lsu2;
