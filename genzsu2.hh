#pragma once

#include<complex>
#include<math.h>
#include<cmath>
#include<cassert>
#include "su2.hh"
#include "parameters.hh"

using Complex = std::complex<double>;

class _Gsu2 {
public:
  const size_t N_c = 2;
  size_t m = 3;
  explicit _Gsu2(){
    j[0] = m;
    j[1] = 0;
    j[2] = 0;
    j[3] = 0;
    s[0] = +1;
    s[1] = +1;
    s[2] = +1;
    s[3] = +1;
  }
  explicit _Gsu2(size_t m) {
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
  explicit _Gsu2(size_t _m,
                 size_t * _j,
                 int * _s){
    for(int i = 0; i < 4; i++) {
      j[i] = _j[i];
      s[i] = _s[i];
    }
    
    assert((m == (j[0] + j[1] + j[2] + j[3])));
    assert(((abs(s[0]) == 1) && (abs(s[1]) == 1) && (abs(s[2]) == 1) && (abs(s[3]) == 1)));
  }
  explicit _Gsu2(size_t m,
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
  /*friend su2 operator*(const _Gsu2 &U1, const _Gsu2 &U2);
  friend su2 operator*(const su2 &U1, const _Gsu2 &U2);
  friend su2 operator*(const _Gsu2 &U1, const su2 &U2);
*/
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
  void operator=(const _Gsu2 &U) {
    for(int i = 0; i < 4; i++) {
      j[i] = U.getj(i);
      s[i] = U.gets(i);
      
    }
    
    assert((m == (j[0] + j[1] + getj(2) + j[3])));
    assert(((abs(s[0]) == 1) && (abs(s[1]) == 1) && (abs(s[2]) == 1) && (abs(s[3]) == 1)));
  }
  
  su2 operator=(const su2 &U){
    return (U);
  }

  su2 set_to_identity() {
    Complex a1(s[0]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[1]*sqrt(static_cast<double>(getj(1))/static_cast<double>(m)));
    Complex b1(s[2]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[3]*sqrt(static_cast<double>(getj(3))/static_cast<double>(m)));
    su2 helpmatrix1(a1, b1);
    return (helpmatrix1 * helpmatrix1.dagger());
  }

  su2 operator*=(const _Gsu2 &U){
  Complex a1(s[0]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[1]*sqrt(static_cast<double>(getj(1))/static_cast<double>(m)));
  Complex b1(s[2]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[3]*sqrt(static_cast<double>(getj(3))/static_cast<double>(m)));
  su2 helpmatrix1(a1, b1);
  Complex a2(U.s[0]*sqrt(static_cast<double>(U.getj(2))/static_cast<double>(m)),
             U.s[1]*sqrt(static_cast<double>(U.getj(1))/static_cast<double>(m)));
  Complex b2(U.s[2]*sqrt(static_cast<double>(U.getj(2))/static_cast<double>(m)),
             U.s[3]*sqrt(static_cast<double>(U.getj(3))/static_cast<double>(m)));
  su2 helpmatrix2(a2, b2);
  return (helpmatrix1 * helpmatrix2);
  }

  /*
  su2 operator+=(const double &U){
    
  Complex a1(s[0]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[1]*sqrt(static_cast<double>(getj(1))/static_cast<double>(m)));
  Complex b1(s[2]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[3]*sqrt(static_cast<double>(getj(3))/static_cast<double>(m)));
  su2 helpmatrix(a1, b1);
  return (helpmatrix + U);
  } */

  su2 operator+=(const su2 &U){
    
  Complex a1(s[0]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[1]*sqrt(static_cast<double>(getj(1))/static_cast<double>(m)));
  Complex b1(s[2]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[3]*sqrt(static_cast<double>(getj(3))/static_cast<double>(m)));
  su2 helpmatrix(a1, b1);
  return (helpmatrix + U);
  }

  /*su2 operator-=(const double &U){
    
  Complex a1(s[0]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[1]*sqrt(static_cast<double>(getj(1))/static_cast<double>(m)));
  Complex b1(s[2]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[3]*sqrt(static_cast<double>(getj(3))/static_cast<double>(m)));
  su2 helpmatrix(a1, b1);
  return (helpmatrix + U);
  } */

  su2 operator-=(const su2 &U){
    
  Complex a1(s[0]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[1]*sqrt(static_cast<double>(getj(1))/static_cast<double>(m)));
  Complex b1(s[2]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[3]*sqrt(static_cast<double>(getj(3))/static_cast<double>(m)));
  su2 helpmatrix(a1, b1);
  return (helpmatrix - U);
  }



  su2 round(size_t n) const {
  Complex a1(s[0]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[1]*sqrt(static_cast<double>(getj(1))/static_cast<double>(m)));
  Complex b1(s[2]*sqrt(static_cast<double>(getj(2))/static_cast<double>(m)),
             s[3]*sqrt(static_cast<double>(getj(3))/static_cast<double>(m)));
  su2 helpmatrix(a1, b1);
  return (helpmatrix.round(n));
  }

  void set(const size_t * const _j, const int * const _s, const size_t _m) {
    for(int i = 0; i < 4; i++) {
      j[i] = _j[i];
      s[i] = _s[i];
    }
    
    assert((m == (j[0] + j[1] + j[2] + j[3])));
    assert(((abs(s[0]) == 1) && (abs(s[1]) == 1) && (abs(s[2]) == 1) && (abs(s[3]) == 1)));
  }
  _Gsu2 dagger() const {
    return(_Gsu2(m, getj(2), getj(1), getj(2), s[0], -s[1], -s[2], -s[3]));
  }
  double retrace() {
    return(2 * s[0] * sqrt(static_cast<double>(getj(2))/static_cast<double>(m)));
  }

  Complex trace(){
    double a = s[0] * sqrt(static_cast<double>(getj(2))/static_cast<double>(m));
    return (Complex(2*a, 0));
  }
  Complex det() {
    return((static_cast<double>(getj(2)) + static_cast<double>(getj(1)) +
            static_cast<double>(getj(2)) + static_cast<double>(getj(3)))/static_cast<double>(m)); 
  }
  void restoreSU() {
    assert((m == (j[0] + j[1] + j[2] + j[3])));
    assert(((abs(s[0]) == 1) && (abs(s[1]) == 1) && (abs(s[2]) == 1) && (abs(s[3]) == 1)));
    
    return;
  }
  double weight() {
    return(1.);
  }

private:
  size_t j[4];
  int s[4];
};

template<> struct accum_type<_Gsu2> {
  typedef _su2 type;
};


su2 operator*(const _Gsu2 &U1, const _Gsu2 &U2) {
  unsigned int m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex a2(U2.gets(0)*sqrt(static_cast<double>(U2.getj(0))/static_cast<double>(m)),
             U2.gets(1)*sqrt(static_cast<double>(U2.getj(1))/static_cast<double>(m)));
  Complex b1(U1.gets(2)*sqrt(static_cast<double>(U1.getj(2))/static_cast<double>(m)),
             U1.gets(3)*sqrt(static_cast<double>(U1.getj(3))/static_cast<double>(m)));
  Complex b2(U2.gets(2)*sqrt(static_cast<double>(U2.getj(2))/static_cast<double>(m)),
             U2.gets(3)*sqrt(static_cast<double>(U2.getj(3))/static_cast<double>(m)));
  su2 res(a1*a2 - b1*std::conj(b2), a1*b2 + b1*std::conj(a2));
  return(res);
}
su2 operator+(const _Gsu2 &U1, const _Gsu2 &U2) {
  unsigned int m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex a2(U2.gets(0)*sqrt(static_cast<double>(U2.getj(0))/static_cast<double>(m)),
             U2.gets(1)*sqrt(static_cast<double>(U2.getj(1))/static_cast<double>(m)));
  Complex b1(U1.gets(2)*sqrt(static_cast<double>(U1.getj(2))/static_cast<double>(m)),
             U1.gets(3)*sqrt(static_cast<double>(U1.getj(3))/static_cast<double>(m)));
  Complex b2(U2.gets(2)*sqrt(static_cast<double>(U2.getj(2))/static_cast<double>(m)),
             U2.gets(3)*sqrt(static_cast<double>(U2.getj(3))/static_cast<double>(m)));
  su2 helpmatrix1(a1, b1);
  su2 helpmatrix2(a2, b2);
  return (helpmatrix1 + helpmatrix2);
}
su2 operator-(const _Gsu2 &U1, const _Gsu2 &U2) {
  unsigned int m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex a2(U2.gets(0)*sqrt(static_cast<double>(U2.getj(0))/static_cast<double>(m)),
             U2.gets(1)*sqrt(static_cast<double>(U2.getj(1))/static_cast<double>(m)));
  Complex b1(U1.gets(2)*sqrt(static_cast<double>(U1.getj(2))/static_cast<double>(m)),
             U1.gets(3)*sqrt(static_cast<double>(U1.getj(3))/static_cast<double>(m)));
  Complex b2(U2.gets(2)*sqrt(static_cast<double>(U2.getj(2))/static_cast<double>(m)),
             U2.gets(3)*sqrt(static_cast<double>(U2.getj(3))/static_cast<double>(m)));
  su2 helpmatrix1(a1, b1);
  su2 helpmatrix2(a2, b2);
  return (helpmatrix1 - helpmatrix2);
}

su2 operator*(const su2 &U1, const _Gsu2 &U2) {
  unsigned int m = U2.getm();
  Complex a2(U2.gets(0)*sqrt(static_cast<double>(U2.getj(0))/static_cast<double>(m)),
             U2.gets(1)*sqrt(static_cast<double>(U2.getj(1))/static_cast<double>(m)));
  Complex b2(U2.gets(2)*sqrt(static_cast<double>(U2.getj(2))/static_cast<double>(m)),
             U2.gets(3)*sqrt(static_cast<double>(U2.getj(3))/static_cast<double>(m)));
  su2 res(U1.geta()*a2 - U1.getb()*std::conj(b2), U1.geta()*b2 + U1.getb()*std::conj(a2));
  return(res);
}
su2 operator+(const _su2 &U1, const _Gsu2 &U2) {
  unsigned int m = U2.getm();
  Complex a2(U2.gets(0)*sqrt(static_cast<double>(U2.getj(0))/static_cast<double>(m)),
             U2.gets(1)*sqrt(static_cast<double>(U2.getj(1))/static_cast<double>(m)));
  Complex b2(U2.gets(2)*sqrt(static_cast<double>(U2.getj(2))/static_cast<double>(m)),
             U2.gets(3)*sqrt(static_cast<double>(U2.getj(3))/static_cast<double>(m)));
  su2 helpmatrix1(a2, b2);
  return (U1 + U2);
}

su2 operator-(const _su2 &U1, const _Gsu2 &U2) {
  unsigned int m = U2.getm();
  Complex a2(U2.gets(0)*sqrt(static_cast<double>(U2.getj(0))/static_cast<double>(m)),
             U2.gets(1)*sqrt(static_cast<double>(U2.getj(1))/static_cast<double>(m)));
  Complex b2(U2.gets(2)*sqrt(static_cast<double>(U2.getj(2))/static_cast<double>(m)),
             U2.gets(3)*sqrt(static_cast<double>(U2.getj(3))/static_cast<double>(m)));
  su2 helpmatrix1(a2, b2);
  return (U1 - U2);
}
su2 operator*(const _Gsu2 &U1, const su2 &U2) {
  unsigned int m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex b1(U1.gets(2)*sqrt(static_cast<double>(U1.getj(2))/static_cast<double>(m)),
             U1.gets(3)*sqrt(static_cast<double>(U1.getj(3))/static_cast<double>(m)));
  su2 res(a1*U2.geta() - b1*std::conj(U2.getb()),  a1*U2.getb() + b1*std::conj(U2.geta()));
  return(res);
}
su2 operator+(const _Gsu2 &U1, const _su2 &U2) {
  unsigned int m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex b1(U1.gets(2)*sqrt(static_cast<double>(U1.getj(2))/static_cast<double>(m)),
             U1.gets(3)*sqrt(static_cast<double>(U1.getj(3))/static_cast<double>(m)));
  su2 helpmatrix2(a1, b1);
  return (U1 + U2);
}

su2 operator-(const _Gsu2 &U1, const _su2 &U2) {
  unsigned int m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex b1(U1.gets(2)*sqrt(static_cast<double>(U1.getj(2))/static_cast<double>(m)),
             U1.gets(3)*sqrt(static_cast<double>(U1.getj(3))/static_cast<double>(m)));
  su2 helpmatrix2(a1, b1);
  return (U1 - U2);
}

su2 operator*(const _Gsu2 &U1, const Complex &U2) {
  size_t m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex b1(U1.gets(2)*sqrt(static_cast<double>(U1.getj(2))/static_cast<double>(m)),
             U1.gets(3)*sqrt(static_cast<double>(U1.getj(3))/static_cast<double>(m)));
  su2 helpmatrix(a1, b1);
  return (helpmatrix * U2); 
}

su2 operator*(const _Gsu2 &U1, const double &U2){
  size_t m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex b1(U1.gets(2)*sqrt(static_cast<double>(U1.getj(2))/static_cast<double>(m)),
             U1.gets(3)*sqrt(static_cast<double>(U1.getj(3))/static_cast<double>(m)));
  su2 helpmatrix(a1, b1);
  return (helpmatrix * U2);
}
su2 operator*(const double &U1, const _Gsu2 &U2) {
  return (U2 *U1);
}
su2 operator*(const Complex &U1, const _Gsu2 &U2){
  return (U2*U1);
}
double retrace(_Gsu2 const &U){
    double a = U.gets(0) * sqrt(static_cast<double>(U.getj(0))/static_cast<double>(U.getm()));

  return (2*a);
}

Complex trace(_Gsu2 const &U){
  double a = U.gets(0) * sqrt(static_cast<double>(U.getj(0))/static_cast<double>(U.getm()));
  return (Complex(2*a, 0));

}


/*inline _Gsu2 accum_to_Group(const Gsu2_accum &x) {
  _Gsu2 U = x;
  U.restoreSU();
  return U;
}


 * @brief accumulation type for SU(2) matrices
 *
 * Incidentally, for SU(2) linear combinations of products its matrices can be still
 * parametrized by 2 numbers only. Therefore we use the same class for SU(2) elements and
 * their accumulations.
 *
 * @tparam
 

template <> struct accum_type<_Gsu2> {
  typedef _Gsu2 type;
};
*/

using Gsu2 = _Gsu2;
