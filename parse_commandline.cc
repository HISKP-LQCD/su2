#include<iostream>
#include<string>
#include<boost/program_options.hpp>

#include"parse_commandline.hh"
#include"print_program_options.hh"

namespace po = boost::program_options;


void add_general_options(po::options_description &desc, general_params &params) {
  desc.add_options()
    ("help,h", "produce this help message")
    ("spatialsizex,X", po::value<size_t>(&params.Lx), "spatial lattice size x > 0")
    ("spatialsizey,Y", po::value<size_t>(&params.Ly), "spatial lattice size y > 0")
    ("spatialsizez,Z", po::value<size_t>(&params.Lz), "spatial lattice size z > 0")
    ("temporalsize,T", po::value<size_t>(&params.Lt), "temporal lattice size > 0")
    ("ndims", po::value<size_t>(&params.ndims), "number of dimensions, 2 <= ndims <= 4")
    ("nsave", po::value<size_t>(&params.N_save)->default_value(1000), "N_save")
    ("n_meas,n", po::value<size_t>(&params.n_meas)->default_value(10), "total number of sweeps")
    ("counter", po::value<size_t>(&params.icounter)->default_value(0), "initial counter for updates")
    ("beta,b", po::value<double>(&params.beta), "beta value")
    ("mass,m", po::value<double>(&params.m0), "bare quark mass")
    ("seed,s", po::value<size_t>(&params.seed)->default_value(13526463), "PRNG seed")
    ("heat", po::value<double>(&params.heat)->default_value(1.), "randomness of the initial config, 1: hot, 0: cold")
    ("restart", "restart from an existing configuration")
    ("configname", po::value< std::string >(&params.configfilename), "configuration filename used in case of restart")
    ("xi", po::value<double>(&params.xi)->default_value(1.0), "xi, characteristic of anisotropy")
    ("anisotropic", po::value<bool>(&params.anisotropic)->default_value(false), "set whether the configurations are anisotropic")
    ;
  return;
}

int parse_commandline(int ac, char * av[], po::options_description &desc, general_params &params) {
  try {
    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    params.restart = false;
    params.acceptreject = true;
    if (vm.count("help")) {
      std::cout << desc << std::endl;
      return 1;
    }
    if (!vm.count("spatialsizex") && !vm.count("help")) {
      std::cerr << "spatial lattice x-size must be given!" << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 1;
    }
    if (!vm.count("spatialsizey") && !vm.count("help")) {
      std::cerr << "spatial lattice y-size must be given!" << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 1;
    }
    if (!vm.count("spatialsizez") && !vm.count("help")) {
      std::cerr << "spatial lattice z-size must be given!" << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 1;
    }
    if (!vm.count("temporalsize") && !vm.count("help")) {
      std::cerr << "temporal lattice size must be given!" << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 1;
    }
    if (!vm.count("beta") && !vm.count("help")) {
      std::cerr << "beta value must be specified!" << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 1;
    }
    if (vm.count("restart")) {
      if(!vm.count("configname")) {
        std::cerr << "for a restart the configuration filename must be specified!" << std::endl;
        std::cout << std::endl << desc << std::endl;
        return 1;
      }
      params.restart = true;
    }
    if (vm.count("no-accept-reject")) {
      params.acceptreject = false;
    }
    if (!vm.count("ndims")) {
      params.ndims = 4;
    }
    if (params.ndims > 4 || params.ndims < 2) {
      std::cerr << "2 <= ndims <= 4!" << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 1;
    }
    if (params.Lx < 1 || params.Ly < 1 || params.Lz < 1 || params.Lt < 1) {
      std::cerr << "All box extends must be > 1!" << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 1;
    }
    if (params.ndims == 2) {
      params.Ly = 1;
      params.Lz = 1;
    }
    if (params.ndims == 3) {
      params.Lz = 1;
    }
    if (params.xi != 1.0 && params.anisotropic == false) {
        std::cerr << "anisotropy parameter xi != 1, but flag anisotropic is set to false. Programm will continue without anisotropy." << std::endl;
    }
    PrintVariableMap(vm);
  }
  catch(std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  catch(...) {
    std::cerr << "Exception of unknown type!\n";
  }
  return(0);
}

int parse_command_line_and_init(int ac, char *av[], general_params& gparams, hmc_u1_params& hparams) {
  po::options_description desc("Allowed options");
  add_general_options(desc, gparams);
  // add HMC specific options
  desc.add_options()("nrev", po::value<size_t>(&hparams.N_rev)->default_value(0),
                     "frequenz of reversibility tests N_rev, 0: not reversibility test")(
    "nsteps", po::value<size_t>(&hparams.n_steps)->default_value(1000), "n_steps")(
    "tau", po::value<double>(&hparams.tau)->default_value(1.), "trajectory length tau")(
    "exponent", po::value<size_t>(&hparams.exponent)->default_value(0),
    "exponent for rounding")("integrator", po::value<size_t>(&hparams.integs)->default_value(0),
                             "itegration scheme to be used: 0=leapfrog, 1=lp_leapfrog, "
                             "2=omf4, 3=lp_omf4, 4=Euler, 5=RUTH, 6=omf2")(
    "no_fermions", po::value<bool>(&hparams.no_fermions)->default_value(0),
    "Bool flag indicating if we're ignoring the fermionic action.")(
    "solver", po::value<std::string>(&hparams.solver)->default_value("CG"),
    "Type of solver: CG, BiCGStab")(
    "tolerace_cg", po::value<double>(&hparams.tolerance_cg)->default_value(1e-10),
    "Tolerance for the solver for the dirac operator")(
    "solver_verbosity", po::value<size_t>(&hparams.solver_verbosity)->default_value(0),
    "Verbosity for the solver for the dirac operator")(
    "seed_pf", po::value<size_t>(&hparams.seed_pf)->default_value(97234719),
    "Seed for the evaluation of the fermion determinant")(
    "outdir", po::value<std::string>(&hparams.outdir)->default_value("."), "Output directory");

  return parse_commandline(ac, av, desc, gparams);
}