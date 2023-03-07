/**
 * @file main-su3.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief main programm running any simulation of this library for the SU(3) theory
 * @version 0.1
 * @date 2022-10-03
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "run_program.hpp"

int main(int argc, char *argv[]) {
  typedef _su3 Group;
  run_program<Group>(argc, argv);

  return (0);
}
