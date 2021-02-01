#include"su2.hh"
#include"ostsu2.hh"
#include<complex>
#include<math.h>

su2 operator*(const Osu2 &U1, const Osu2 &U2) {
  const double J1 = U1.getJ(), J2 = U2.getJ();
  Complex a1(U1.s[0]*static_cast<double>(U1.j[0])/J1,
             U1.s[1]*static_cast<double>(U1.j[1])/J1);
  Complex a2(U2.s[0]*static_cast<double>(U2.j[0])/J2,
             U2.s[1]*static_cast<double>(U2.j[1])/J2);
  Complex b1(U1.s[2]*static_cast<double>(U1.j[2])/J1,
             U1.s[3]*static_cast<double>(U1.j[3])/J1);
  Complex b2(U2.s[2]*static_cast<double>(U2.j[2])/J2,
             U2.s[3]*static_cast<double>(U2.j[3])/J2);
  su2 res(a1*a2 - b1*std::conj(b2), a1*b2 + b1*std::conj(a2));
  return(res);
}

su2 operator*(const su2 &U1, const Osu2 &U2) {
  const double J2 = U2.getJ();
  Complex a2(U2.s[0]*static_cast<double>(U2.j[0])/J2,
             U2.s[1]*static_cast<double>(U2.j[1])/J2);
  Complex b2(U2.s[2]*static_cast<double>(U2.j[2])/J2,
             U2.s[3]*static_cast<double>(U2.j[3])/J2);
  su2 res(U1.geta()*a2 - U1.getb()*std::conj(b2), U1.geta()*b2 + U1.getb()*std::conj(a2));
  return(res);
}

su2 operator*(const Osu2 &U1, const su2 &U2) {
  const double J1 = U1.getJ();
  Complex a1(U1.s[0]*static_cast<double>(U1.j[0])/J1,
             U1.s[1]*static_cast<double>(U1.j[1])/J1);
  Complex b1(U1.s[2]*static_cast<double>(U1.j[2])/J1,
             U1.s[3]*static_cast<double>(U1.j[3])/J1);
  su2 res(a1*U2.geta() - b1*std::conj(U2.getb()),  a1*U2.getb() + b1*std::conj(U2.geta()));
  return(res);
}
