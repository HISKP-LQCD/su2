// parse_input_file.cc

#include <iostream>
#include <string>

#include "parse_input_file.hh"

namespace YAML_parsing {

  /**
   * @brief saving value from YAML node
   * Saving in 'x' the value specified in the YAML node under the key string 'name'.
   * If the key doesn't exist, nothing is done
   * @tparam T cast type of the parameter
   * @param x address of the value
   * @param nd YAML node
   * @param name string name of the parameter
   */
  template <class T> void read(T &x, const YAML::Node &nd, const std::string &name) {    
    try {
      x = nd[name].as<T>();
    } catch(...) {
      std::cerr << "Error: check \""<< name << "\" in your YAML input file. ";
      std::cerr << boost::typeindex::type_id<T>() << " type was expected. \n";
      abort();
    }
    return;
  }

  // read() and output on std::cout
  template <class T> void read_verbose(T &x, const YAML::Node &nd, const std::string &name) {    
    read(x, nd, name);
    std::cout << "## " << name << "=" << x << "\n";
    return;
  }

  /* 
  Reading optional argument: same as read(), but if the key doesn't exist, nothing is done
  This function should be used with parameters that have a default argument.
  */
  template <class T>
  void read_optional(T &x, const YAML::Node &nd, const std::string &name) {
    if (nd[name]) {
      read<T>(x, nd, name);
    }
    return;
  }

  // read_optional() and output on std::cout
  template <class T>
  void read_optional_verbose(T &x, const YAML::Node &nd, const std::string &name) {
    if (nd[name]) {
      read<T>(x, nd, name);
    }
    else{
      std::cout << "## " << name << "=" << x << " (default)\n";
    }
    return;
  }

} // namespace YAML_parsing

int parse_input_file(const std::string &file,
                     gp::general &gparams,
                     gp::hmc_u1 &hmc_params) {

  const YAML::Node nd = YAML::LoadFile(file);

  namespace Yp = YAML_parsing;
  std::cout << "## Parsing input file: " << file << "\n";

  // general parameters
  std::cout << nd["geometry"]["X"] << "\n";
  Yp::read_verbose<size_t>(gparams.Lx, nd["geometry"], "X");
  Yp::read_verbose<size_t>(gparams.Ly, nd["geometry"], "Y");
  Yp::read_verbose<size_t>(gparams.Lz, nd["geometry"], "Z");
  Yp::read_verbose<size_t>(gparams.Lt, nd["geometry"], "T");
  Yp::read_verbose<size_t>(gparams.ndims, nd["geometry"], "ndims");
  
  Yp::read_verbose<double>(gparams.beta, nd["action"], "beta");
  Yp::read_verbose<double>(gparams.m0, nd["action"], "mass");
  
  Yp::read_verbose<size_t>(gparams.N_save, nd["hmc_trajectories"], "nsave");
  Yp::read_verbose<size_t>(gparams.N_meas, nd["hmc_trajectories"], "nmeas");
  Yp::read_verbose<size_t>(gparams.icounter, nd["hmc_trajectories"], "counter");
  Yp::read_verbose<size_t>(gparams.seed, nd["hmc_trajectories"], "seed");
  Yp::read_verbose<double>(gparams.heat, nd["hmc_trajectories"], "heat");
  Yp::read_verbose<std::string>(gparams.configfilename, nd["hmc_trajectories"], "configname");

  if (gparams.ndims > 4 || gparams.ndims < 2) {
    std::cerr << "2 <= ndims <= 4!" << std::endl;
    return 1;
  }
  if (gparams.Lx < 1 || gparams.Ly < 1 || gparams.Lz < 1 || gparams.Lt < 1) {
    std::cerr << "All box extends must be > 1!" << std::endl;
    return 1;
  }
  if (gparams.ndims == 2) {// flattening 'y' and 'z' directions
    gparams.Ly = 1;
    gparams.Lz = 1;
  }
  if (gparams.ndims == 3) {// flattening 'z' direction
    gparams.Lz = 1;
  }

  // hmc-u1 parameters
  Yp::read_optional_verbose<size_t>(hmc_params.N_rev, nd["hmc_trajectories"], "N_rev");
  Yp::read_optional_verbose<size_t>(hmc_params.n_steps, nd["hmc_trajectories"], "n_steps");
  Yp::read_optional_verbose<double>(hmc_params.tau, nd["hmc_trajectories"], "tau");
  Yp::read_optional_verbose<size_t>(hmc_params.exponent, nd["hmc_trajectories"], "exponent");
  Yp::read_optional_verbose<size_t>(hmc_params.integs, nd["hmc_trajectories"], "integs");
  Yp::read_optional_verbose<bool>(hmc_params.no_fermions, nd["hmc_trajectories"], "no_fermions");
  Yp::read_optional_verbose<std::string>(hmc_params.solver, nd["hmc_trajectories"], "solver");
  Yp::read_optional_verbose<double>(hmc_params.tolerance_cg, nd["hmc_trajectories"], "tolerance_cg");
  Yp::read_optional_verbose<size_t>(hmc_params.solver_verbosity, nd["hmc_trajectories"], "solver_verbosity");
  Yp::read_optional_verbose<size_t>(hmc_params.seed_pf, nd["hmc_trajectories"], "seed_pf");
  Yp::read_optional_verbose<std::string>(hmc_params.outdir, nd["hmc_trajectories"], "outdir");

  return 0;
}
