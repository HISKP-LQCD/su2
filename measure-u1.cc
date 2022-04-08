#include"su2.hh"
#include"u1.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"sweep.hh"
#include"wilsonloop.hh"
#include"md_update.hh"
#include"monomial.hh"
#include"gradient_flow.hh"
#include"energy_density.hh"
#include"parse_input_file.hh"
#include"version.hh"

#include<iostream>
#include<iomanip>
#include<sstream>
#include<vector>
#include<random>
#include<boost/program_options.hpp>

namespace po = boost::program_options;

#include <boost/filesystem.hpp>

int main(int ac, char* av[]) {
  std::cout << "## Measuring Tool for U(1) gauge theory" << std::endl;
  std::cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017)" << std::endl;
  std::cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << std::endl;

  namespace gp = global_parameters;
  gp::physics pparams; // physics parameters
  gp::measure_u1 mparams; // measure parameters

  std::string input_file; // yaml input file path
  po::options_description desc("Allowed options");
  desc.add_options()
  ("help,h", "produce this help message")
  ("file,f", po::value<std::string>(&input_file)->default_value("NONE"), "yaml input file");

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 0;
  }

  namespace in_meas = input_file_parsing::u1::measure;
  int err = in_meas::parse_input_file(input_file, pparams, mparams);
  if (err > 0) {
    return err;
  }

  boost::filesystem::create_directories(boost::filesystem::absolute(mparams.confdir));

  gaugeconfig<_u1> U(pparams.Lx, pparams.Ly, pparams.Lz, pparams.Lt, pparams.ndims, pparams.beta);

    std::stringstream ss_basename;
    ss_basename << mparams.conf_basename << ".";
    ss_basename << pparams.Lx << "." << pparams.Ly << "." << pparams.Lz << "."
                << pparams.Lt;
    ss_basename << ".b" << std::fixed << std::setprecision(mparams.beta_str_width)
                << pparams.beta;
  for(size_t i = mparams.icounter; i < mparams.n_meas*mparams.nstep+mparams.icounter; i+=mparams.nstep) {
    std::ostringstream os;
    os << mparams.confdir << ss_basename.str() << "." << i << std::ends;
    int ierrU =  U.load(os.str());
    if(ierrU == 1){ // cannot load gauge config
      continue;
    }
    
    double plaquette = gauge_energy(U);
    double density = 0., Q=0.;
    energy_density(U, density, Q);
    std::cout << "## Initial Plaquette: " << plaquette/U.getVolume()/double(U.getNc())/6. << std::endl; 
    std::cout << "## Initial Energy density: " << density << std::endl;
    
    random_gauge_trafo(U, mparams.seed);
    plaquette = gauge_energy(U);
    energy_density(U, density, Q);
    std::cout << "## Plaquette after rnd trafo: " << std::scientific << std::setw(15) << plaquette/U.getVolume()/double(U.getNc())/6. << std::endl; 
    std::cout << "## Energy density: " << density << std::endl;
        
    if(mparams.Wloop) {
      std::ostringstream os;
      os << mparams.confdir + "/wilsonloop.";
      auto prevw = os.width(6);
      auto prevf = os.fill('0');
      os << i;
      os.width(prevw);
      os.fill(prevf);
      os << ".dat" << std::ends;
      compute_all_loops(U, os.str());
    }
    if(mparams.gradient) {
      std::ostringstream os;
      os << mparams.confdir + "/gradient_flow.";
      auto prevw = os.width(6);
      auto prevf = os.fill('0');
      os << i;
      os.width(prevw);
      os.fill(prevf);
      gradient_flow(U, os.str(), mparams.tmax);
    }
  }

  return(0);
}
