/**
 * @file hmc-u1.cc
 * @author Carsten Urbach (urbach@hiskp.uni-bonn.de)
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief Hybrid Monte Carlo for a U(1) theory
 * @version 0.1
 * @date 2022-05-11
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "hmc-u1.hpp"

int main(int ac, char *av[]) {

  u1::hmc_algo h1;
  h1.run(ac, av);

  return (0);
}
