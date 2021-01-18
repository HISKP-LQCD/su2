#include"su2.hh"
#include"genzsu2.hh"
#include<complex>
#include<math.h>

su2 operator*(const Gsu2 &U1, const Gsu2 &U2) {
  unsigned int m = U1.m;
  Complex a1(U1.s[0]*sqrt(static_cast<double>(U1.j[0])/static_cast<double>(m)),
             U1.s[1]*sqrt(static_cast<double>(U1.j[1])/static_cast<double>(m)));
  Complex a2(U2.s[0]*sqrt(static_cast<double>(U2.j[0])/static_cast<double>(m)),
             U2.s[1]*sqrt(static_cast<double>(U2.j[1])/static_cast<double>(m)));
  Complex b1(U1.s[2]*sqrt(static_cast<double>(U1.j[2])/static_cast<double>(m)),
             U1.s[3]*sqrt(static_cast<double>(U1.j[3])/static_cast<double>(m)));
  Complex b2(U2.s[2]*sqrt(static_cast<double>(U2.j[2])/static_cast<double>(m)),
             U2.s[3]*sqrt(static_cast<double>(U2.j[3])/static_cast<double>(m)));
  su2 res(a1*a2 - b1*std::conj(b2), a1*b2 + b1*std::conj(a2));
  return(res);
}

su2 operator*(const su2 &U1, const Gsu2 &U2) {
  unsigned int m = U2.m;
  Complex a2(U2.s[0]*sqrt(static_cast<double>(U2.j[0])/static_cast<double>(m)),
             U2.s[1]*sqrt(static_cast<double>(U2.j[1])/static_cast<double>(m)));
  Complex b2(U2.s[2]*sqrt(static_cast<double>(U2.j[2])/static_cast<double>(m)),
             U2.s[3]*sqrt(static_cast<double>(U2.j[3])/static_cast<double>(m)));
  su2 res(U1.geta()*a2 - U1.getb()*std::conj(b2), U1.geta()*b2 + U1.getb()*std::conj(a2));
  return(res);
}

su2 operator*(const Gsu2 &U1, const su2 &U2) {
  unsigned int m = U1.m;
  Complex a1(U1.s[0]*sqrt(static_cast<double>(U1.j[0])/static_cast<double>(m)),
             U1.s[1]*sqrt(static_cast<double>(U1.j[1])/static_cast<double>(m)));
  Complex b1(U1.s[2]*sqrt(static_cast<double>(U1.j[2])/static_cast<double>(m)),
             U1.s[3]*sqrt(static_cast<double>(U1.j[3])/static_cast<double>(m)));
  su2 res(a1*U2.geta() - b1*std::conj(U2.getb()),  a1*U2.getb() + b1*std::conj(U2.geta()));
  return(res);
}
