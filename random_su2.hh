#include<random>
#include"su2.hh"

template<class URNG> void random_su2(su2 &U, URNG &engine, const double delta) {
  std::uniform_real_distribution<double> dist1(-1., 1.);
  std::uniform_real_distribution<double> dist2(0., 2*3.1415);
  std::uniform_real_distribution<double> dist3(0., delta*2*3.1415);
  const double alpha = dist3(engine);
  const double u = dist1(engine);
  const double theta = dist2(engine);
  
  const double r = sqrt(1-u*u);
  const double salpha = sin(alpha);

  U = su2(Complex(cos(alpha), salpha*u), salpha*r*Complex(sin(theta), cos(theta)));
  return;
}
