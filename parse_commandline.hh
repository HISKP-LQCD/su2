#pragma once
#include<string>
#include<boost/program_options.hpp>
namespace po = boost::program_options;

class general_params {
public:
  size_t Ls, Lt;             
  size_t N_meas;             
  size_t N_save;             
  size_t icounter;
  size_t seed;                  
  double beta;               
  double heat;               
  bool restart;              
  std::string configfilename;
};

void add_general_options(po::options_description &desc, general_params &params);
int parse_commandline(int ac, char * av[], po::options_description &desc, general_params &params);
