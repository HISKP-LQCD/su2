// main-u1.cc
/**
 * @file main-u1.cc
 * @author Carsten Urbach (urbach@hiskp.uni-bonn.de)
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief Metropolis Algorithm for U(1) gauge theory
 * @version 0.1
 * @date 2022-05-09
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "metropolis-u1.hpp"

int main(int ac, char *av[]) {

  u1::metropolis m1;
  m1.run(ac, av);

  return (0);
}
