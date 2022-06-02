#pragma once

#include <boost/program_options.hpp>
#include <string>

// #include "parameters.hh"

namespace po = boost::program_options;

// namespace gp = global_parameters;
// typedef general_params general_params;

struct general_params {
public:
  size_t Lx; // spatial lattice size x > 0
  size_t Ly; // spatial lattice size y > 0
  size_t Lz; // spatial lattice size z > 0
  size_t Lt; // temporal lattice size > 0
  size_t ndims = 4; // number of dimensions, 2 <= ndims <= 4
  double beta; // beta value
  double m0; // bare quark mass
  double xi = 1.0; // anisotropy parameter

  size_t N_save = 1000; // N_save
  size_t n_meas = 10; // total number of sweeps
  size_t icounter = 0; // initial counter for updates
  size_t seed = 13526463; // PRNG seed
  double heat = 1.0; // randomness of the initial config, 1: hot, 0: cold
  bool restart = false; // restart from an existing configuration
  bool acceptreject = true; // use accept-reject
  bool anisotropic = false; // use anisotropic lattice
  std::string configfilename = ""; // configuration filename used in case of restart
};

struct hmc_u1_params {
  size_t N_rev = 0; // frequency of reversibility tests N_rev, 0: no reversibility test
  size_t n_steps = 1000; // n_steps
  double tau = 1; // trajectory length tau
  size_t exponent = 0; // exponent for rounding
  size_t integs = 0; // itegration scheme to be used: 0=leapfrog, 1=lp_leapfrog, 2=omf4,
                     // 3=lp_omf4, 4=Euler, 5=RUTH, 6=omf2
  bool no_fermions = 0; // Bool flag indicating if we're ignoring the fermionic action.
  std::string solver = "CG"; // Type of solver: CG, BiCGStab
  double tolerance_cg = 1e-10; // Tolerance for the solver for the dirac operator
  size_t solver_verbosity = 0; // Verbosity for the solver for the dirac operator
  size_t seed_pf = 97234719; // Seed for the evaluation of the fermion determinant
  std::string outdir = "."; // Output directory
};

//~ struct measure_u1_params {
    //~ size_t nmeas; //number of configurations that should be measured
    //~ size_t nstep = 1; //measure each nstep config
    //~ bool Wloop = false; //measure Wilson-loops
    //~ bool gradient_flow= false; //measure Grandient flow
    //~ double tmax = 9.99; //tmax for gradient flow
    //~ std::string confdir = "."; //directory where configs are stored
    //~ bool potential = false; //measure potential: the loops W(x, t, y=z=0) and W(x, y, t=z=0) are measured with a maximum size of lattice extent * sizeloops, and written to separate files. Only available for ndims=3,4
    //~ bool potentialsmall = false; //The loops W(x, t, y) are measured up to x, y=min(4, lattice extent), t <= Lt * sizeloops and saved to one file. Only available for ndim=3
    //~ bool append = false; //are measurements for potential appended to an existing file, or should it be overwritten?
    //~ double sizeloops = 0.5; //Wilson-Loops are measured up to this fraction of the lattice extent
    //~ size_t n_apesmear = 0; //number of APE smearings done on the lattice before measurement. 
    //~ //APE-smearing is done before measuring the potential and small potential, it does not affect the gradient flow and Wilson-loops
    //~ double alpha = 1.0; //parameter alpha for APE smearings
    //~ bool spacialsmear = false; //should smearing be done only for spacial links?
    //~ // missing tau, exponent, integs
//~ };

void add_general_options(po::options_description &desc, general_params &params);
int parse_commandline(int ac,
                      char *av[],
                      po::options_description &desc,
                      general_params &params);

// parsing the command line and initializing the values in gparams and hmc_params
int parse_command_line_and_init(int ac,
                                char *av[],
                                general_params &gparams,
                                hmc_u1_params &hparams);
