#pragma once

#include"su2.hh"
#include"genzsu2.hh"
#include"ostsu2.hh"
#include<random>
#include<iostream>

constexpr double pi() { return std::atan(1.)*4.; }

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

template<class URNG> void random_su2(Gsu2 &U, URNG &engine, 
                                     const size_t m = 10,
                                     const double delta = 0.) {

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
}

template<class URNG> void random_su2(Osu2 &U, URNG &engine, 
                                     const size_t m = 10,
                                     const double delta = 0.) {


  size_t lower = static_cast<size_t>(delta * m);
  size_t j[4];
  std::uniform_int_distribution<int> uni1(lower, m);
  j[0] = uni1(engine);
  std::uniform_int_distribution<int> uni2(0, m-j[0]);
  j[1] = uni2(engine);
  std::uniform_int_distribution<int> uni3(0, m-j[0]-j[1]);
  j[2] = uni3(engine);
  j[3] = m - j[0] - j[1] - j[2];

  std::uniform_int_distribution<int> uni4(0, 1);
  int s[4] = {1,1,1,1};
  int ilower = static_cast<int>(delta * 4);
  for(int i = ilower; i < 4; i++) {
    s[i] = 2*uni4(engine) - 1;
  }
  U = Osu2(m, j, s);
  return;
}
