#pragma once

#include<complex>
#include<math.h>
#include<cmath>
#include<cassert>
#include"su2.hh"

using Complex = std::complex<double>;

class _Gsu2 {
public:
  explicit _Gsu2() : m(0) {
    j[0] = m;
    j[1] = 0;
    j[2] = 0;
    j[3] = 0;
    s[0] = +1;
    s[1] = +1;
    s[2] = +1;
    s[3] = +1;
  }
  explicit _Gsu2(unsigned int m) : m(m) {
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
  explicit _Gsu2(unsigned int _m,
                 unsigned int * _j,
                 int * _s){
    for(int i = 0; i < 4; i++) {
      j[i] = _j[i];
      s[i] = _s[i];
    }
    m=_m;
    assert((m == (j[0] + j[1] + j[2] + j[3])));
    assert(((abs(s[0]) == 1) && (abs(s[1]) == 1) && (abs(s[2]) == 1) && (abs(s[3]) == 1)));
  }
  explicit _Gsu2(unsigned int m,
                 unsigned int j1, unsigned int j2, unsigned int j3,
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
  friend su2 operator*(const _Gsu2 &U1, const _Gsu2 &U2);
  friend su2 operator*(const su2 &U1, const _Gsu2 &U2);
  friend su2 operator*(const _Gsu2 &U1, const su2 &U2);

  unsigned int getj(const unsigned int i) const {
    return(j[i]);
  }
  int gets(const unsigned int i) const {
    return(s[i]);
  }
  unsigned int getm() const {
    return(m);
  }
  void operator=(const _Gsu2 &U) {
    for(int i = 0; i < 4; i++) {
      j[i] = U.getj(i);
      s[i] = U.gets(i);
      m = U.getm();
    }
    assert((m == (j[0] + j[1] + j[2] + j[3])));
    assert(((abs(s[0]) == 1) && (abs(s[1]) == 1) && (abs(s[2]) == 1) && (abs(s[3]) == 1)));
  }
  void set(const unsigned int * const _j, const int * const _s, const unsigned int _m) {
    for(int i = 0; i < 4; i++) {
      j[i] = _j[i];
      s[i] = _s[i];
    }
    m=_m;
    assert((m == (j[0] + j[1] + j[2] + j[3])));
    assert(((abs(s[0]) == 1) && (abs(s[1]) == 1) && (abs(s[2]) == 1) && (abs(s[3]) == 1)));
  }
  _Gsu2 dagger() const {
    return(_Gsu2(m, j[0], j[1], j[2], s[0], -s[1], -s[2], -s[3]));
  }
  double trace() {
    return(2 * s[0] * sqrt(static_cast<double>(j[0])/static_cast<double>(m)));
  }
  Complex det() {
    return((static_cast<double>(j[0]) + static_cast<double>(j[1]) +
            static_cast<double>(j[2]) + static_cast<double>(j[3]))/static_cast<double>(m)); 
  }
  void restoreSU() {
    return;
  }

private:
  unsigned int m;
  unsigned int j[4];
  int s[4];
};

su2 operator*(const _Gsu2 &U1, const _Gsu2 &U2);
su2 operator*(const su2 &U1, const _Gsu2 &U2);
su2 operator*(const _Gsu2 &U1, const su2 &U2);

using Gsu2 = _Gsu2;
