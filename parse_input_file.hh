#pragma once

#include <boost/type_index.hpp>

#include "parameters.hh"
#include "yaml.h"

namespace gp = global_parameters;


/**
 * @brief input file parsing
 * parsing the input file, assigning the values to the attributes of gparams and hmc_params
 * return
 * @param file input file path
 * @param gparams global parameter structure
 * @param hmc_params hmc parameter structure 
 * @return int error value
 */
int parse_input_file(const std::string &file,
                      gp::general &gparams,
                      gp::hmc_u1 &hmc_params);
