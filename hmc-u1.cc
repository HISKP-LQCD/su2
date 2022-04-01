#include"su2.hh"
#include"u1.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"energy_density.hh"
#include"random_gauge_trafo.hh"
#include"md_update.hh"
#include"monomial.hh"
#include"integrator.hh"
#include"parse_commandline.hh"
#include"parse_input_file.hh"
#include"version.hh"

#include<iostream>
#include<fstream>
#include<iomanip>
#include<sstream>
#include<random>
#include<boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "detDDdag_monomial.hh"

namespace po = boost::program_options;
using std::cout;
using std::endl;

int main(int ac, char *av[]) {

  cout << "## HMC Algorithm for U(1) gauge theory" << endl;
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017, 2021)" << endl;
  cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << endl;

  namespace gp = global_parameters;
  gp::general gparams; // general parameters
  gp::hmc_u1 hmc_params; // hmc parameters

  std::string input_file; // yaml input file path
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce this help message")(
    "file,f", po::value<std::string>(&input_file)->default_value("NONE"),
    "yaml input file");
  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << desc << "\n";
    return 0;
  }

  int err = parse_input_file(input_file, gparams, hmc_params);
  if (err > 0) {
    return 1;
  }

  boost::filesystem::create_directories(boost::filesystem::absolute(hmc_params.outdir));
  
  gaugeconfig<_u1> U(gparams.Lx, gparams.Ly, gparams.Lz, gparams.Lt, gparams.ndims,
                     gparams.beta);
  if (gparams.restart) {
    err = U.load(gparams.configfilename);
    if (err != 0) {
      return err;
    }
  } else {
    hotstart(U, gparams.seed, gparams.heat);
  }

  double plaquette = gauge_energy(U);
  double fac = 2. / U.getndims() / (U.getndims() - 1);
  const double normalisation = fac / U.getVolume() / double(U.getNc());
  cout << "## Initital Plaquette: " << plaquette * normalisation << endl;

  random_gauge_trafo(U, 654321);
  plaquette = gauge_energy(U);
  cout << "## Plaquette after rnd trafo: " << plaquette * normalisation << endl;

  // Molecular Dynamics parameters
  md_params mdparams(hmc_params.n_steps, hmc_params.tau);

  // generate list of monomials
  gaugemonomial<double, _u1> gm(0);
  kineticmonomial<double, _u1> km(0);
  detDDdag_monomial<double, _u1> detDDdag(0, gparams.m0, hmc_params.solver,
                                          hmc_params.tolerance_cg, hmc_params.seed_pf,
                                          hmc_params.solver_verbosity);

  km.setmdpassive();

  std::list<monomial<double, _u1> *> monomial_list;
  monomial_list.push_back(&gm);
  monomial_list.push_back(&km);

  if (!hmc_params.no_fermions) { // including S_F (fermionic) in the action
    monomial_list.push_back(&detDDdag);
  }

  integrator<double, _u1> *md_integ = set_integrator<double, _u1>(hmc_params.integs, hmc_params.exponent);

  std::ofstream os;
  if (gparams.icounter == 0)
    os.open(hmc_params.outdir+"/output.hmc.data", std::ios::out);
  else
    os.open(hmc_params.outdir+"/output.hmc.data", std::ios::app);

  std::cout << "## Normalization factor: A = 2/(d*(d-1)*N_lat*N_c) = " << std::scientific
            << std::setw(18) << std::setprecision(15) << normalisation << "\n";
  std::cout << "## Acceptance rate parcentage: rho = rate/(i+1)\n";

  std::stringstream ss_head; // header: column names in the output
  ss_head<<"i"<<" "<<"getaccept"<<" "<< "E*A"<<" "<<"dH"<<" "<<"rho"<<" "<<"ddH"<<" "<<"Q"<<"\n";

  std::string ss_head_str = ss_head.str();
  std::cout << ss_head_str;
  os << ss_head_str;

  double rate = 0.;
  for (size_t i = gparams.icounter; i < gparams.N_meas + gparams.icounter; i++) {
    mdparams.disablerevtest();
    if (i > 0 && hmc_params.N_rev != 0 && (i) % hmc_params.N_rev == 0) {
      mdparams.enablerevtest();
    }
    // PRNG engine
    std::mt19937 engine(gparams.seed + i);
    // perform the MD update
    md_update(U, engine, mdparams, monomial_list, *md_integ);

    double energy = gauge_energy(U);
    double E = 0., Q = 0.;
    energy_density(U, E, Q);
    rate += mdparams.getaccept();

    cout << i << " " << mdparams.getaccept() << " " << std::scientific << std::setw(18)
         << std::setprecision(15) << energy * normalisation << " " << std::setw(15)
         << mdparams.getdeltaH() << " " << std::setw(15)
         << rate / static_cast<double>(i + 1) << " ";
    if (mdparams.getrevtest()) {
      cout << mdparams.getdeltadeltaH();
    } else
      cout << "NA";
    cout << " " << Q << endl;

    os << i << " " << mdparams.getaccept() << " " << std::scientific << std::setw(18)
       << std::setprecision(15) << energy * normalisation << " " << std::setw(15)
       << mdparams.getdeltaH() << " " << std::setw(15)
       << rate / static_cast<double>(i + 1) << " ";
    if (mdparams.getrevtest()) {
      os << mdparams.getdeltadeltaH();
    } else
      os << "NA";
    os << " " << Q << endl;

    if (i > 0 && (i % gparams.N_save) == 0) {// saving U after each N_save trajectories
      std::ostringstream oss;
      oss << "config_u1." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "."
          << gparams.Lt << ".b" << gparams.beta << "." << i << std::ends;
      U.save(hmc_params.outdir+"/"+oss.str());
    }
  }
  cout << "## Acceptance rate: " << rate / static_cast<double>(gparams.N_meas) << endl;

  std::ostringstream oss;
  oss << "config_u1." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "."
      << gparams.Lt << ".b" << U.getBeta() << ".final" << std::ends;
  U.save(hmc_params.outdir+"/"+oss.str());
  return (0);
}
