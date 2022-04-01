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

  // reading optional argument: same as read(), but if the key doesn't exist, nothing is done
  template <class T>
  void read_optional(T &x, const YAML::Node &nd, const std::string &name) {
    if (nd[name]) {
      read<T>(x, nd, name);
    }
    return;
  }

} // namespace YAML_parsing

int parse_input_file(const std::string &file,
                     gp::general &gparams,
                     gp::hmc_u1 &hmc_params) {
  // try {
  //   po::variables_map vm;
  //   po::store(po::parse_command_line(ac, av, desc), vm);
  //   po::notify(vm);

  //   params.restart = false;
  //   params.acceptreject = true;
  //   if (vm.count("help")) {
  //     std::cout << desc << std::endl;
  //     return 1;
  //   }
  //   if (!vm.count("spatialsizex") && !vm.count("help")) {
  //     std::cerr << "spatial lattice x-size must be given!" << std::endl;
  //     std::cout << std::endl << desc << std::endl;
  //     return 1;
  //   }
  //   if (!vm.count("spatialsizey") && !vm.count("help")) {
  //     std::cerr << "spatial lattice y-size must be given!" << std::endl;
  //     std::cout << std::endl << desc << std::endl;
  //     return 1;
  //   }
  //   if (!vm.count("spatialsizez") && !vm.count("help")) {
  //     std::cerr << "spatial lattice z-size must be given!" << std::endl;
  //     std::cout << std::endl << desc << std::endl;
  //     return 1;
  //   }
  //   if (!vm.count("temporalsize") && !vm.count("help")) {
  //     std::cerr << "temporal lattice size must be given!" << std::endl;
  //     std::cout << std::endl << desc << std::endl;
  //     return 1;
  //   }
  //   if (!vm.count("beta") && !vm.count("help")) {
  //     std::cerr << "beta value must be specified!" << std::endl;
  //     std::cout << std::endl << desc << std::endl;
  //     return 1;
  //   }
  //   if (vm.count("restart")) {
  //     if(!vm.count("ndname")) {
  //       std::cerr << "for a restart the nduration filename must be specified!" <<
  //       std::endl; std::cout << std::endl << desc << std::endl; return 1;
  //     }
  //     params.restart = true;
  //   }
  //   if (vm.count("no-accept-reject")) {
  //     params.acceptreject = false;
  //   }
  //   if (!vm.count("ndims")) {
  //     params.ndims = 4;
  //   }
  //   if (params.ndims > 4 || params.ndims < 2) {
  //     std::cerr << "2 <= ndims <= 4!" << std::endl;
  //     std::cout << std::endl << desc << std::endl;
  //     return 1;
  //   }
  //   if (params.Lx < 1 || params.Ly < 1 || params.Lz < 1 || params.Lt < 1) {
  //     std::cerr << "All box extends must be > 1!" << std::endl;
  //     std::cout << std::endl << desc << std::endl;
  //     return 1;
  //   }
  //   if (params.ndims == 2) {
  //     params.Ly = 1;
  //     params.Lz = 1;
  //   }
  //   if (params.ndims == 3) {
  //     params.Lz = 1;
  //   }
  //   PrintVariableMap(vm);
  // }
  // catch(std::exception& e) {
  //   std::cerr << "error: " << e.what() << "\n";
  //   return 1;
  // }
  // catch(...) {
  //   std::cerr << "Exception of unknown type!\n";
  // }

  const YAML::Node nd = YAML::LoadFile(file);

  namespace Yp = YAML_parsing;

  // general parameters
  Yp::read<size_t>(gparams.Lx, nd, "X");
  Yp::read<size_t>(gparams.Ly, nd, "Y");
  Yp::read<size_t>(gparams.Lz, nd, "Z");
  Yp::read<size_t>(gparams.Lt, nd, "T");
  Yp::read<size_t>(gparams.ndims, nd, "ndims");
  Yp::read<size_t>(gparams.N_save, nd, "nsave");
  Yp::read<size_t>(gparams.N_meas, nd, "nmeas");
  Yp::read<size_t>(gparams.icounter, nd, "counter");
  Yp::read<double>(gparams.beta, nd, "beta");
  Yp::read<double>(gparams.m0, nd, "mass");
  Yp::read<size_t>(gparams.seed, nd, "seed");
  Yp::read<double>(gparams.heat, nd, "heat");
  Yp::read<std::string>(gparams.configfilename, nd, "configname");

  // hmc-u1 parameters
  Yp::read_optional<size_t>(hmc_params.N_rev, nd, "N_rev");
  Yp::read_optional<size_t>(hmc_params.n_steps, nd, "n_steps");
  Yp::read_optional<double>(hmc_params.tau, nd, "tau");
  Yp::read_optional<size_t>(hmc_params.exponent, nd, "exponent");
  Yp::read_optional<size_t>(hmc_params.integs, nd, "integs");
  Yp::read_optional<bool>(hmc_params.no_fermions, nd, "no_fermions");
  Yp::read_optional<std::string>(hmc_params.solver, nd, "solver");
  Yp::read_optional<double>(hmc_params.tolerance_cg, nd, "tolerance_cg");
  Yp::read_optional<size_t>(hmc_params.solver_verbosity, nd, "solver_verbosity");
  Yp::read_optional<size_t>(hmc_params.seed_pf, nd, "seed_pf");
  Yp::read_optional<std::string>(hmc_params.outdir, nd, "outdir");

  return (0);
}
