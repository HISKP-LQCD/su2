// parse_input_file.hh
/**
 * @brief parsing input file routines
 * This file contains the following namespaces:
 * - YAML_parsing : parsing parameters from YAML nodes
 * - input_file_parsing : parsing sets of parameters from files
 */

#pragma once

#include <boost/type_index.hpp>

#include "parameters.hh"
#include "yaml.h"

/**
 * @brief namespace of yaml parsing functions
 * This namespace contains functions that read parameters specified in the input file from
 * which YAML nodes are generated
 */
namespace YAML_parsing {

  /**
   * @brief saving value from YAML node
   * Saving in 'x' the value specified in the YAML node under the key string 'name'.
   * If the key doesn't exist, nothing is done
   * @tparam T cast type of the parameter
   * @param x address of the value
   * @param nd YAML node
   * @param name string name of the parameter
   */
  template <class T> void read(T &x, const YAML::Node &nd, const std::string &name);

  // read() and output on std::cout
  template <class T> void read_verb(T &x, const YAML::Node &nd, const std::string &name);

  /*
  Reading optional argument: same as read(), but if the key doesn't exist, nothing is done
  This function should be used with parameters that have a default argument.
  */
  template <class T>
  void read_optional(T &x, const YAML::Node &nd, const std::string &name);

  // read_optional() and output on std::cout
  template <class T>
  void read_opt_verb(T &x, const YAML::Node &nd, const std::string &name);

} // namespace YAML_parsing

/**
 * @brief parsing from input files
 * This namespace defines functions that, starting from a given input file, initialize the
 * attributes of structures needed for the various programs. The structures are passed as
 * reference to the functions. Their structure is:
 * @param file input file path
 * @param pparams physics parameter structure
 * @param xxx_params xxx program-specific parameter structure (e.g. xxx=hmc_params)
 * @return int error value
 */
namespace input_file_parsing {

  namespace gp = global_parameters;

  /**
   * @brief check geometry parameters
   * Check if the lattice spacetime extensions are sensible and consistent with the number
   * of dimensions (d). It is checked 2<=d. For 2<d<4, spatial dimensions of 'pparams' are
   * flattened according to the following priority order: Lz, Ly, Lx
   * @param file
   * @param pparams
   * @return int
   */
  int validate_geometry(gp::physics &pparams);

  namespace u1 {

  /**
   * @brief parsing geometrical parameters
   * Parsing the spacetime lattice extensions and number of dimensions from the input file.
   * If invalid, the program aborts
   * @param pparams physics parameter structure
   */
    void parse_geometry(const YAML::Node &nd, gp::physics &pparams);

    namespace hmc {
      int parse_input_file(const std::string &file,
                           gp::physics &pparams,
                           gp::hmc_u1 &hmc_params);
    }

    namespace measure {
      int parse_input_file(const std::string &file,
                           gp::physics &pparams,
                           gp::measure_u1 &mparams);

    } // namespace measure

  } // namespace u1

} // namespace input_file_parsing
