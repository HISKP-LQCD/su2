/**
 * @file main-su2.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief main programm running any simulation of this library for the SU(2) theory
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

#include "hmc.hpp"
#include "measure.hpp"
#include "metropolis.hpp"

typedef _su2 Group;

int main(int argc, char *argv[]) {
  std::string input_file;
  parse_command_line(argc, argv, input_file);
  std::cout << "## Parsing input file: " << input_file << "\n";

  running_program rp; // running program
  YAML::Node nd = YAML::Clone(get_cleaned_input_file(rp, input_file));

  bool &do_hmc = rp.do_hmc;
  bool &do_metropolis = rp.do_metropolis;
  bool &do_omeas = rp.do_omeas;

if (do_hmc ^ do_metropolis) { // one of the 2 algorithms
    if (do_hmc) {
      hmc_algo<Group> h;
      h.run(nd);
    } else if (do_metropolis) {
      metropolis_algo<Group> mpl;
      mpl.run(nd);
    }
  } else if (do_omeas) { // offline measurements
    measure_algo<Group> ms;
    ms.run(nd);
  } else { // program does nothing
    std::cerr << "ERROR: Program ineffective.\n";
    return 1;
  }

  return (0);
}
