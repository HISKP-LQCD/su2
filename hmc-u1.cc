#include"su2.hh"
#include"u1.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"energy_density.hh"
#include"random_gauge_trafo.hh"
#include"md_update.hh"
#include"monomial.hh"
#include"integrator.hh"
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
using std::endl;

int main(int ac, char *av[]) {

  std::cout << "## HMC Algorithm for U(1) gauge theory" << std::endl;
  std::cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017, 2021)" << std::endl;
  std::cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << std::endl;

  namespace gp = global_parameters;
  gp::physics pparams; // physics parameters
  gp::hmc_u1 hparams; // hmc parameters

  std::string input_file; // yaml input file path
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce this help message")(
    "file,f", po::value<std::string>(&input_file)->default_value("NONE"),
    "yaml input file");
  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 0;
  }

  namespace in_hmc = input_file_parsing::u1::hmc;
  int err = in_hmc::parse_input_file(input_file, pparams, hparams);
  if (err > 0) {
    return 1;
  }

  boost::filesystem::create_directories(boost::filesystem::absolute(hparams.outdir));
  
  gaugeconfig<_u1> U(pparams.Lx, pparams.Ly, pparams.Lz, pparams.Lt, pparams.ndims,
                     pparams.beta);
  if (hparams.restart) {
    std::cout << "restart " <<hparams.restart << "\n";
    err = U.load(hparams.configfilename);
    if (err != 0) {
      return err;
    }
  } else {
    std::cout << "hotstart "<<hparams.seed<<" "<< hparams.heat<<"\n";
    const double heat_val = (hparams.heat==true)? 1.0: 0.0;
    hotstart(U, hparams.seed, heat_val);
  }

  double plaquette = gauge_energy(U);
  double fac = 2. / U.getndims() / (U.getndims() - 1);
  const double normalisation = fac / U.getVolume() / double(U.getNc());
  std::cout << "## Initital Plaquette: " << plaquette * normalisation << std::endl;

  random_gauge_trafo(U, 654321);
  plaquette = gauge_energy(U);
  std::cout << "## Plaquette after rnd trafo: " << plaquette * normalisation << std::endl;

  // Molecular Dynamics parameters
  md_params mdparams(hparams.n_steps, hparams.tau);

  // generate list of monomials
  std::list<monomial<double, _u1> *> monomial_list;
  gaugemonomial<double, _u1> gm(0);

  kineticmonomial<double, _u1> km(0);
  km.setmdpassive();
  monomial_list.push_back(&km);

  detDDdag_monomial<double, _u1> detDDdag(0, pparams.m0, hparams.solver,
                                            hparams.tolerance_cg, hparams.seed_pf,
                                            hparams.solver_verbosity);


  if (pparams.include_gauge) {
    monomial_list.push_back(&gm);
  }

  if (pparams.include_staggered_fermions) { // including S_F (fermionic) in the action
    monomial_list.push_back(&detDDdag);
  }


  // setting up the integrator
  integrator<double, _u1> *md_integ = set_integrator<double, _u1>(hparams.integrator, hparams.exponent);

  std::ofstream os;
  if (hparams.icounter == 0)
    os.open(hparams.outdir+"/output.hmc.data", std::ios::out);
  else
    os.open(hparams.outdir+"/output.hmc.data", std::ios::app);

  std::cout << "## Normalization factor: A = 2/(d*(d-1)*N_lat*N_c) = " << std::scientific
            << std::setw(18) << std::setprecision(15) << normalisation << "\n";
  std::cout << "## Acceptance rate parcentage: rho = rate/(i+1)\n";

  std::stringstream ss_head; // header: column names in the output
  ss_head<<"i"<<" "<<"getaccept"<<" "<< "E*A"<<" "<<"dH"<<" "<<"rho"<<" "<<"ddH"<<" "<<"Q"<<"\n";

  std::string ss_head_str = ss_head.str();
  std::cout << ss_head_str;
  os << ss_head_str;

  double rate = 0.;

  const std::string conf_basename = hparams.conf_basename;
  std::stringstream ss_basename;
  ss_basename << hparams.conf_basename << ".";
  ss_basename << pparams.Lx << "." << pparams.Ly << "." << pparams.Lz << "."
              << pparams.Lt;
  ss_basename << ".b" << std::fixed << std::setprecision(hparams.beta_str_width)
              << pparams.beta;

  for (size_t i = hparams.icounter; i < hparams.n_meas + hparams.icounter; i++) {
    mdparams.disablerevtest();
    if (i > 0 && hparams.N_rev != 0 && (i) % hparams.N_rev == 0) {
      mdparams.enablerevtest();
    }
    // PRNG engine
    std::mt19937 engine(hparams.seed + i);
    // perform the MD update
    md_update(U, engine, mdparams, monomial_list, *md_integ);

    double energy = gauge_energy(U);
    double E = 0., Q = 0.;
    energy_density(U, E, Q);
    rate += mdparams.getaccept();

    std::cout << i << " " << mdparams.getaccept() << " " << std::scientific << std::setw(18)
         << std::setprecision(15) << energy * normalisation << " " << std::setw(15)
         << mdparams.getdeltaH() << " " << std::setw(15)
         << rate / static_cast<double>(i + 1) << " ";
    if (mdparams.getrevtest()) {
      std::cout << mdparams.getdeltadeltaH();
    } else
      std::cout << "NA";
    std::cout << " " << Q << std::endl;

    os << i << " " << mdparams.getaccept() << " " << std::scientific << std::setw(18)
       << std::setprecision(15) << energy * normalisation << " " << std::setw(15)
       << mdparams.getdeltaH() << " " << std::setw(15)
       << rate / static_cast<double>(i + 1) << " ";
    if (mdparams.getrevtest()) {
      os << mdparams.getdeltadeltaH();
    } else
      os << "NA";
    os << " " << Q << std::endl;

    if (i > 0 && (i % hparams.N_save) == 0) { // saving U after each N_save trajectories
      std::ostringstream oss_i;
      oss_i << ss_basename.str() << "." << i << std::ends;
      U.save(hparams.outdir + "/" + oss_i.str());
    }
  }
  std::cout << "## Acceptance rate: " << rate / static_cast<double>(hparams.n_meas)
            << std::endl;

  std::ostringstream oss;
  oss << ss_basename.str() << ".final" << std::ends;
  U.save(hparams.outdir + "/" + oss.str());

  return (0);
}
