/**
 * @file main-u1.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief main programm running any simulation of this library for the U(1) theory
 * @version 0.1
 * @date 2022-10-03
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <ctime>
#include <iostream>
#include <sstream>

#include "boost/lexical_cast.hpp"

#include "parse_input_file.hh"

#include "hmc-u1.hpp"
#include "measure-u1.hpp"
#include "metropolis-u1.hpp"

int main(int argc, char *argv[]) {
  std::string input_file;
  u1::parse_command_line(argc, argv, input_file);
  std::cout << "## Parsing input file: " << input_file << "\n";

  u1::running_program rp; // running program
  YAML::Node nd = YAML::Clone(u1::get_cleaned_input_file(rp, input_file));

  typedef _u1 Group;
  std::cout << rp.gg << "-ciaoooo\n";
  if (rp.gg == "u1") {
    typedef _u1 Group;
  }
  else if (rp.gg == "su2") {
    typedef _su2 Group;
  } else {
    spacetime_lattice::fatal_error(
      "Invalid gauge group specified in the input file: " + rp.gg, __func__);
  }

  bool &do_hmc = rp.do_hmc;
  bool &do_metropolis = rp.do_metropolis;
  bool &do_omeas = rp.do_omeas;

if (do_hmc ^ do_metropolis) { // one of the 2 algorithms
    if (do_hmc) {
      u1::hmc_algo<Group> h;
      h.run(nd);
    } else if (do_metropolis) {
      u1::metropolis_algo<Group> mpl;
      mpl.run(nd);
    }
  } else if (do_omeas) { // offline measurements
    u1::measure_algo<Group> ms;
    ms.run(nd);
  } else { // program does nothing
    std::cerr << "ERROR: Program ineffective.\n";
    return 1;
  }

  return (0);
}
