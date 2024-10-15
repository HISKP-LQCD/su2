#include"su2.hh"
#include"genzsu2.hh"
#include<complex>
#include<math.h>


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

su2 operator-(const Gsu2 &U1, const Gsu2 &U2) {
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

su2 operator*(const Gsu2 &U1, const Gsu2 &U2) {
  unsigned int m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex a2(U2.gets(0)*sqrt(static_cast<double>(U2.getj(0))/static_cast<double>(m)),
             U2.gets(1)*sqrt(static_cast<double>(U2.getj(1))/static_cast<double>(m)));
  Complex b1(U1.s[2]*sqrt(static_cast<double>(U1.j[2])/static_cast<double>(m)),
             U1.s[3]*sqrt(static_cast<double>(U1.j[3])/static_cast<double>(m)));
  Complex b2(U2.s[2]*sqrt(static_cast<double>(U2.j[2])/static_cast<double>(m)),
             U2.s[3]*sqrt(static_cast<double>(U2.j[3])/static_cast<double>(m)));
  su2 res(a1*a2 - b1*std::conj(b2), a1*b2 + b1*std::conj(a2));
  return(res);
}

su2 operator*(const su2 &U1, const Gsu2 &U2) {
  unsigned int m = U2.getm();
  Complex a2(U2.gets(0)*sqrt(static_cast<double>(U2.getj(0))/static_cast<double>(m)),
             U2.gets(1)*sqrt(static_cast<double>(U2.getj(1))/static_cast<double>(m)));
  Complex b2(U2.s[2]*sqrt(static_cast<double>(U2.j[2])/static_cast<double>(m)),
             U2.s[3]*sqrt(static_cast<double>(U2.j[3])/static_cast<double>(m)));
  su2 res(U1.geta()*a2 - U1.getb()*std::conj(b2), U1.geta()*b2 + U1.getb()*std::conj(a2));
  return(res);
}

su2 operator*(const Gsu2 &U1, const su2 &U2) {
  unsigned int m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex b1(U1.s[2]*sqrt(static_cast<double>(U1.j[2])/static_cast<double>(m)),
             U1.s[3]*sqrt(static_cast<double>(U1.j[3])/static_cast<double>(m)));
  su2 res(a1*U2.geta() - b1*std::conj(U2.getb()),  a1*U2.getb() + b1*std::conj(U2.geta()));
  return(res);
}

su2 operator+(const su2 &U1, const Gsu2 &U2) {
  unsigned int m = U2.getm();
  Complex a2(U2.gets(0)*sqrt(static_cast<double>(U2.getj(0))/static_cast<double>(m)),
             U2.gets(1)*sqrt(static_cast<double>(U2.getj(1))/static_cast<double>(m)));
  Complex b2(U2.s[2]*sqrt(static_cast<double>(U2.j[2])/static_cast<double>(m)),
             U2.s[3]*sqrt(static_cast<double>(U2.j[3])/static_cast<double>(m)));
  su2 helpmatrix1(a2, b2);
  return (U1 + U2);
}

su2 operator+(const Gsu2 &U1, const su2 &U2) {
  unsigned int m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex b1(U1.s[2]*sqrt(static_cast<double>(U1.j[2])/static_cast<double>(m)),
             U1.s[3]*sqrt(static_cast<double>(U1.j[3])/static_cast<double>(m)));
  su2 helpmatrix2(a1, b1);
  return (U1 + U2);
}

su2 operator-(const su2 &U1, const Gsu2 &U2) {
  unsigned int m = U2.getm();
  Complex a2(U2.gets(0)*sqrt(static_cast<double>(U2.getj(0))/static_cast<double>(m)),
             U2.gets(1)*sqrt(static_cast<double>(U2.getj(1))/static_cast<double>(m)));
  Complex b2(U2.s[2]*sqrt(static_cast<double>(U2.j[2])/static_cast<double>(m)),
             U2.s[3]*sqrt(static_cast<double>(U2.j[3])/static_cast<double>(m)));
  su2 helpmatrix1(a2, b2);
  return (U1 - U2);
}

su2 operator-(const Gsu2 &U1, const su2 &U2) {
  unsigned int m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex b1(U1.s[2]*sqrt(static_cast<double>(U1.j[2])/static_cast<double>(m)),
             U1.s[3]*sqrt(static_cast<double>(U1.j[3])/static_cast<double>(m)));
  su2 helpmatrix2(a1, b1);
  return (U1 - U2);
}




su2 operator*(const _Gsu2 &U1, const Complex &U2) {
  size_t m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex b1(U1.s[2]*sqrt(static_cast<double>(U1.j[2])/static_cast<double>(m)),
             U1.s[3]*sqrt(static_cast<double>(U1.j[3])/static_cast<double>(m)));
  su2 helpmatrix(a1, b1);
  return (helpmatrix * U2); 
}

su2 operator*(const _Gsu2 &U1, const double &U2){
  size_t m = U1.getm();
  Complex a1(U1.gets(0)*sqrt(static_cast<double>(U1.getj(0))/static_cast<double>(m)),
             U1.gets(1)*sqrt(static_cast<double>(U1.getj(1))/static_cast<double>(m)));
  Complex b1(U1.s[2]*sqrt(static_cast<double>(U1.j[2])/static_cast<double>(m)),
             U1.s[3]*sqrt(static_cast<double>(U1.j[3])/static_cast<double>(m)));
  su2 helpmatrix(a1, b1);
  return (helpmatrix * U2);
}

su2 operator*(const Complex &U1, const _Gsu2 &U2){
  return (U2*U1)
}