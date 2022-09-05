// measure-u1.cc
/**
 * @file measure-u1.cc
 * @author Carsten Urbach (urbach@hiskp.uni-bonn.de)
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief offline measurements of observables over previously generate gauge
 * configurations
 * @version 0.1
 * @date 2022-09-05
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "measure-u1.hpp"

int main(int ac, char *av[]) {

  u1::measure_algo m2;
  m2.run(ac, av);

  return (0);
}
