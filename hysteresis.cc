#include"su2.hh"
#include"linearsu2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"energy_density.hh"
#include"sweep.hh"
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

  size_t N_hit = 10, dN_hit = 10;
  size_t delta = 2;
  size_t m = 100;

  cout << "## Hysteresis loop using the Metropolis Algorithm for SU(2) gauge theory" << endl;
  cout << "## based on a special partitioning for SU(2)" << endl;
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017-2021)" << endl;
  cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << endl << endl;  

  po::options_description desc("Allowed options");
  add_general_options(desc, gparams);

  // add Metropolis specific options
  desc.add_options()
    ("nhit", po::value<size_t>(&N_hit)->default_value(10), "N_hit")
    ("Genzm,m", po::value<size_t>(&m)->default_value(100), "Genz m")
    ("delta,d", po::value<size_t>(&delta)->default_value(2), "delta")
    ("dnhit", po::value<size_t>(&dN_hit)->default_value(10), "dN_hit")
    ;

  int err = parse_commandline(ac, av, desc, gparams);
  if(err > 0) {
    return err;
  }
  gaugeconfig<Lsu2> U(gparams.Lx, gparams.Ly, gparams.Lz, gparams.Lt, gparams.ndims, gparams.beta);
  for(size_t i = 0; i < U.getSize(); i++) U[i].setm(m);
  hotstart(U, gparams.seed, 0);

  double plaquette = gauge_energy(U);
  double fac = 2./U.getndims()/(U.getndims()-1);
  const double normalisation = fac/U.getVolume()/double(U.getNc());
  cout << "Initital Plaquette: " << plaquette*normalisation << endl; 

  std::ofstream os;
  os.open("output.hysteresis.hot.data", std::ios::out);

  double rate = 0.;
  for(double beta = 0.1; beta <= gparams.beta; beta += 0.1) {
    U.setBeta(beta);
    for(size_t i = gparams.icounter; i < gparams.N_meas + gparams.icounter; i++) {
      std::mt19937 engine(gparams.seed+i);
      rate += sweep(U, engine, m, delta, N_hit, dN_hit, beta);
    }
    double energy = gauge_energy(U);
    double density, Q;
    energy_density(U, density, Q);
    cout << beta << " " << std::scientific << std::setw(5) << std::setprecision(5) << energy*normalisation << endl;
    os << beta << " " << std::scientific << std::setw(5) << std::setprecision(5) << energy*normalisation << endl;
  }
  os.close();
  
  hotstart(U, gparams.seed, 1);

  os.open("output.hysteresis.cold.data", std::ios::out);

  rate = 0.;
  for(double beta = gparams.beta; beta > 0.05; beta -= 0.1) {
    U.setBeta(beta);
    for(size_t i = gparams.icounter; i < gparams.N_meas + gparams.icounter; i++) {
      std::mt19937 engine(gparams.seed+i);
      rate += sweep(U, engine, m, delta, N_hit, dN_hit, beta);
    }
    double energy = gauge_energy(U);
    double density, Q;
    energy_density(U, density, Q);
    cout << beta << " " << std::scientific << std::setw(5) << std::setprecision(5) << energy*normalisation << endl;
    os << beta << " " << std::scientific << std::setw(5) << std::setprecision(5) << energy*normalisation << endl;
  }
  os.close();
  
  return(0);
}
