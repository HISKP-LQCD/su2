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

  u1::running_program rp; // running program
  YAML::Node nd = YAML::Clone(u1::get_cleaned_input_file(rp, input_file));
  std::string tif = u1::get_exported_node_timestamp(nd, input_file); // timestaped input file

  bool &do_hmc = rp.do_hmc;
  bool &do_metropolis = rp.do_metropolis;
  bool &do_omeas = rp.do_omeas;

  if (do_hmc && do_metropolis) { // both options are incompatible
    std::cerr << "ERROR: Can't run simultaneously hmc and metropolis algorithms.\n";
    std::cerr << "Check your input file: " << input_file << "\n";
  } else if (do_hmc ^ do_metropolis) { // one of the 2 algorithms
    if (do_hmc) {
      u1::hmc_algo h;
      h.run(tif);
    } else if (do_metropolis) {
      u1::metropolis_algo mpl;
      mpl.run(tif);
    }
  } else if (do_omeas) { // offline measurements
    u1::measure_algo ms;
    ms.run(tif);
  } else { // program does nothing
    std::cerr << "ERROR: Program ineffective.\n";
    return 1;
  }

  return (0);
}
