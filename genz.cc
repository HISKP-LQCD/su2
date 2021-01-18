#include"su2.hh"
#include"genzsu2.hh"
// #include"gaugeconfig.hh"
// #include"gauge_energy.hh"
// #include"random_gauge_trafo.hh"
// #include"sweep.hh"
#include"parse_commandline.hh"
#include"version.hh"

#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<vector>
#include<random>
#include<boost/program_options.hpp>

using std::vector;
using std::cout;
using std::endl;
namespace po = boost::program_options;

int main(int ac, char* av[]) {
  general_params gparams;

  size_t N_hit = 10;
  double delta = 0.1;

  cout << "## Metropolis Algorithm for SU(2) gauge theory" << endl;
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017,2020,2021)" << endl;
  cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << endl << endl;  

  po::options_description desc("Allowed options");
  add_general_options(desc, gparams);

  // add Metropolis specific options
  desc.add_options()
    ("nhit", po::value<size_t>(&N_hit)->default_value(10), "N_hit")
    ("delta,d", po::value<double>(&delta), "delta")
    ;

  int err = parse_commandline(ac, av, desc, gparams);
  if(err > 0) {
    return err;
  }

  unsigned int j[4] = {3, 2, 4, 1};
  int s[4] = {1, 1, 1, 1};
  
  Gsu2 U1(10, j, s), U2(10);
  su2 U = U1*U2;
  su2 UU = U*U1*U2*U1;
  cout << U1.det() << " " << U1.trace() << endl;
  cout << U2.det() << " " << U2.trace() << endl;
  cout << UU.det() << " " << UU.trace() << endl;
  
  return(0);
}
