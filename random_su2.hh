#include"su2.hh"
#include<random>

constexpr double pi() { return std::atan(1)*4; }

template<class URNG> void random_su2(su2 &U, URNG &engine, 
                                     const double delta = 1.) {

  std::uniform_real_distribution<double> dist1(-1., 1.);
  std::uniform_real_distribution<double> dist2(0., 2*pi());
  std::uniform_real_distribution<double> dist3(0., delta*2*pi());
  const double alpha = dist3(engine);
  const double u = dist1(engine);
  const double theta = dist2(engine);
  
  const double r = sqrt(1-u*u);
  const double salpha = sin(alpha);

  U = su2(Complex(cos(alpha), salpha*u), 
          salpha*r*Complex(sin(theta), cos(theta)));
  return;
}
