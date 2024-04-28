/**
 * @file main-u1.hpp
 * @author Simone Romiti (simone.romiti.1994@gmail.com)
 * @brief main program running any simulation of this library for the U(1) theory
 * @version 0.1
 * @date 2022-10-03
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "run_program.hpp"

int main(int argc, char *argv[]) {
  typedef _u1 Group;
  run_program<Group>(argc, argv);

  return (0);
}
