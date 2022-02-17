#pragma once
#include<string>
#include<boost/program_options.hpp>
namespace po = boost::program_options;

class general_params {
public:
  size_t Lx, Ly, Lz, Lt, ndims;             
  size_t N_meas;             
  size_t N_save;             
  size_t icounter;
  size_t seed; 
  double beta;               
  double heat;               
  bool restart;              
  bool acceptreject;
  std::string configfilename;
  double m0; // quark bare mass (in lattice units)
};

void add_general_options(po::options_description &desc, general_params &params);
int parse_commandline(int ac, char * av[], po::options_description &desc, general_params &params);
