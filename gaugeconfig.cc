#include<random>
#include<vector>
#include<cmath>
#include"su2.hh"
#include"gaugeconfig.hh"


gaugeconfig coldstart(size_t Ls, size_t Lt) {

  gaugeconfig config(Ls, Lt);
  for(int i = 0; i < config.getSize(); i++) {
    config[i] = su2(1., 0.);
  }
  return(config);
}
gaugeconfig hotstart(size_t Ls, size_t Lt, 
                     const int seed, const double _delta) {

  gaugeconfig config(Ls, Lt);
  double delta = _delta;
  if(delta < 0.) delta = 0;
  if(delta > 1.) delta = 1.;
  std::mt19937 engine(seed);
  std::uniform_real_distribution<double> dist1(-1., 1.);
  std::uniform_real_distribution<double> dist2(0., 2*3.1415);
  std::uniform_real_distribution<double> dist3(0., delta*3.1415/2.);

  for(int i = 0; i < config.getSize(); i++) {
    const double alpha = dist3(engine);
    const double u = dist1(engine);
    const double theta = dist2(engine);

    const double r = sqrt(1-u*u);
    const double salpha = sin(alpha);

    config[i] = su2(Complex(cos(alpha), salpha*u), salpha*r*Complex(sin(theta), cos(theta)));
  }
  return(config);
}
