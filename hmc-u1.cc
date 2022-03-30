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
  general_params gparams;
  size_t N_rev;
  size_t n_steps;
  size_t exponent;
  double tau;
  size_t integs;
  bool no_fermions;
  std::string solver;
  double tolerance_cg;
  size_t solver_verbosity;
  size_t seed_pf;
  std::string outdir;

  cout << "## HMC Algorithm for U(1) gauge theory" << endl;
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017, 2021)" << endl;
  cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << endl
       << endl;

  po::options_description desc("Allowed options");
  add_general_options(desc, gparams);
  // add HMC specific options
  desc.add_options()("nrev", po::value<size_t>(&N_rev)->default_value(0),
                     "frequenz of reversibility tests N_rev, 0: not reversibility test")(
    "nsteps", po::value<size_t>(&n_steps)->default_value(1000), "n_steps")(
    "tau", po::value<double>(&tau)->default_value(1.), "trajectory length tau")(
    "exponent", po::value<size_t>(&exponent)->default_value(0),
    "exponent for rounding")("integrator", po::value<size_t>(&integs)->default_value(0),
                             "itegration scheme to be used: 0=leapfrog, 1=lp_leapfrog, "
                             "2=omf4, 3=lp_omf4, 4=Euler, 5=RUTH, 6=omf2")(
    "no_fermions", po::value<bool>(&no_fermions)->default_value(0),
    "Bool flag indicating if we're ignoring the fermionic action.")(
    "solver", po::value<std::string>(&solver)->default_value("CG"),
    "Type of solver: CG, BiCGStab")(
    "tolerace_cg", po::value<double>(&tolerance_cg)->default_value(1e-10),
    "Tolerance for the solver for the dirac operator")(
    "solver_verbosity", po::value<size_t>(&solver_verbosity)->default_value(0),
    "Verbosity for the solver for the dirac operator")(
    "seed_pf", po::value<size_t>(&seed_pf)->default_value(97234719),
    "Seed for the evaluation of the fermion determinant")(
    "outdir", po::value<std::string>(&outdir)->default_value("."),
    "Output directory");

  int err = parse_commandline(ac, av, desc, gparams);
  if (err > 0) {
    return err;
  }

  boost::filesystem::create_directories(outdir);
  
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
  md_params mdparams(n_steps, tau);

  // generate list of monomials
  gaugemonomial<double, _u1> gm(0);
  kineticmonomial<double, _u1> km(0);
  detDDdag_monomial<double, _u1> detDDdag(0, gparams.m0, solver, tolerance_cg, seed_pf,
                                             solver_verbosity);

  km.setmdpassive();

  std::list<monomial<double, _u1> *> monomial_list;
  monomial_list.push_back(&gm);
  monomial_list.push_back(&km);

  if (!no_fermions) { // not neglecting S_F
    monomial_list.push_back(&detDDdag);
  }

  integrator<double, _u1> *md_integ = set_integrator<double, _u1>(integs, exponent);

  std::ofstream os;
  if (gparams.icounter == 0)
    os.open(outdir+"/output.hmc.data", std::ios::out);
  else
    os.open(outdir+"/output.hmc.data", std::ios::app);

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
    if (i > 0 && N_rev != 0 && (i) % N_rev == 0) {
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
    os << " " << Q << " " << endl;

    if (i > 0 && (i % gparams.N_save) == 0) {// saving U after each N_save trajectories
      std::ostringstream oss;
      oss << "config_u1." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "."
          << gparams.Lt << ".b" << gparams.beta << "." << i << std::ends;
      U.save(outdir+"/"+oss.str());
    }
  }
  cout << "## Acceptance rate: " << rate / static_cast<double>(gparams.N_meas) << endl;

  std::ostringstream oss;
  oss << "config_u1." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "."
      << gparams.Lt << ".b" << U.getBeta() << ".final" << std::ends;
  U.save(outdir+"/"+oss.str());
  return (0);
}
