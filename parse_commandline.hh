#pragma once

#include<string>
#include<boost/program_options.hpp>

#include "parameters.hh"

namespace po = boost::program_options;

namespace gp = global_parameters;
typedef gp::general general_params;

void add_general_options(po::options_description &desc, gp::general &params);
int parse_commandline(int ac, char * av[], po::options_description &desc, gp::general &params);


// parsing the command line and initializing the values in gparams and hmc_params
int parse_command_line_and_init(int ac, char *av[], gp::general& gparams, gp::hmc_u1& hmc_params);

