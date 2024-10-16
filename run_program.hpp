/**
 * @file program.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief Running program function. Template argument specified gauge group
 * @version 0.1
 * @date 2022-10-16
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#
#include <ctime>
#include <iostream>
#include <sstream>


#include "boost/lexical_cast.hpp"

#include "parse_input_file.hh"

//#ifndef Genz
#ifndef parti
#include "hmc.hpp"
#include "heatbath_overrelaxation.hpp"
#endif
//#endif

#include "measure.hpp"
#include "metropolis.hpp"

/**
 * @brief program function for hmc, metropolis and measure
 * This function takes the `int main()` function parameters and runs the MCMC simulation.
 * Depending on the input file passed through `argv` it decides whether to
 * - do the MCMC with metropolis or hmc
 * - do the offline/online measurements
 *
 * @tparam Group gauge group: _u1 or _su2
 * @param argc
 * @param argv
 */
template <class Group> void run_program(int argc, char *argv[]) {
  
  std::string input_file;
  parse_command_line(argc, argv, input_file);
  std::cout << "## Parsing input file: " << input_file << "\n";

  running_program rp; // running program
  YAML::Node nd = YAML::Clone(get_cleaned_input_file(rp, input_file));

  bool &do_hmc = rp.do_hmc;
  bool &do_metropolis = rp.do_metropolis;
  bool &do_omeas = rp.do_omeas;
  bool &do_heatbath_overrelaxation = rp.do_heatbath_overrelaxation;

  
  if (do_metropolis) {
    metropolis_algo<Group> mpl;
    mpl.run(nd);
  } 
  //#ifndef Genz
  #ifndef parti
  else if (do_hmc ) {
    hmc_algo<Group> h;
    h.run(nd);
  }
  
   else if (do_heatbath_overrelaxation ) {
    heatbath_overrelaxation_algo<Group> hb_or;
    hb_or.run(nd);
  }
  //#endif
  #endif
   else if (do_omeas) { // offline measurements
    measure_algo<Group> ms;
    ms.run(nd);
  } else { // program does nothing
    std::cerr << "ERROR: Program ineffective.\n";
    exit(1);
  }

  return;
}
