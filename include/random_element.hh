#pragma once

#include "su2.hh"
#include "su3.hh"
#include "u1.hh"
//#include "genzsu2.hh"
#include "partitionings.hh"

#include <random>

constexpr double pi() {
  return std::atan(1) * 4;
}

template <class URNG, class T>
void random_element(T &U, URNG &engine, const double delta = 1.);

template <class URNG> void random_element(_u1 &U, URNG &engine, const double delta = 1.) {
  std::uniform_real_distribution<double> dist(-pi() * delta, pi() * delta);

  U = _u1(dist(engine));
  return;
}

template <class URNG>
void random_element(_su2 &U, URNG &engine, const double delta = 1.) {
  std::uniform_real_distribution<double> dist1(-1., 1.);
  std::uniform_real_distribution<double> dist2(0., 2 * pi());
  std::uniform_real_distribution<double> dist3(0., delta * 2 * pi());
  const double alpha = dist3(engine);
  const double u = dist1(engine);
  const double theta = dist2(engine);

  const double r = sqrt(1 - u * u);
  const double salpha = sin(alpha);

  U = _su2(Complex(cos(alpha), salpha * u), salpha * r * Complex(sin(theta), cos(theta)));
  return;
}

/**
 * @brief Initizalizes the configuration of a random partitioning element with distance smaller than 
 * delta to the identity
 * 
 * @tparam URNG 
 * @param U group element that is replaced with a random element 
 * @param engine engine for random number generation
 * @param delta maximal distance to identity
 */
template <class URNG>
void random_element(partitioning &U, URNG &engine, const double delta = 1.){
  std::uniform_int_distribution<int> dist1(0, partitioning::distance_to_identity.size() -1 );
  // The if statement protects against infinite searches for delta smaller than the smallest distance
  if (delta <=  partitioning::min_distance){
    U.set_to_identity();
  }
  else{
  size_t r = dist1(engine);
  while (_partitioning::distance_to_identity[r] > delta ){
    r = dist1(engine);
  } 
  U = _partitioning(r);}
  return;
}

/**
 * @brief initialized the configuration to some random SU(3) matrix
 *
 * In the (u,v) representation of U (see eq. 4.26 of Gattringer&Lang), there's no unique
 * way of defining "v", because the orthogonal space to a vector "u" is 2-dimensional.
 * Here we adopt our custom "u_2"-convention, according to which: 
 * v = (u_2, -u_1-u_3, u_2)^{\dagger} 
 * At the end we normalize the vectors using the restoreSU() method.
 *
 * @tparam URNG
 * @param U
 * @param engine
 * @param delta
 */
template <class URNG>
void random_element(_su3 &U, URNG &engine, const double delta = 1.) {

  std::uniform_real_distribution<double> dist(-delta, delta);
  std::array<Complex, 3> u;
  for (size_t i = 0; i < U.N_c; i++) {
    u[i] = dist(engine);
  }
  std::array<Complex, 3> v = {std::conj(u[1]), -std::conj(u[0]) - std::conj(u[2]),
                              std::conj(u[1])};
  U = _su3(u, v);
  U.restoreSU();

  return;
}
/* 
template<class URNG> void random_element(Gsu2 &U, URNG &engine, 
                                         const double delta = 1.) {
  const size_t m = U.getm();
  size_t lower = static_cast<int>(delta * m);
  size_t j[4];
  std::uniform_int_distribution<int> uni1(lower, m);
  j[0] = uni1(engine);
  std::uniform_int_distribution<int> uni2(0, m-j[0]);
  j[1] = uni2(engine);
  std::uniform_int_distribution<int> uni3(0, m-j[0]-j[1]);
  j[2] = uni3(engine);
  j[3] = m - j[0] - j[1] - j[2];

  std::uniform_int_distribution<int> uni4(0, 1);
  int s[4];
  for(int i = 0; i < 4; i++) {
    s[i] = 2*uni4(engine) - 1;
  }
  U = Gsu2(m, j, s);
  return;
}*/