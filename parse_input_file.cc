// parse_input_file.cc

#include <iostream>
#include <string>

#include "parse_input_file.hh"

namespace YAML_parsing {

  template <class T> void read(T &x, const YAML::Node &nd, const std::string &name) {
    try {
      x = nd[name].as<T>();
    } catch (...) {
      std::cerr << "Error: check \"" << name << "\" in your YAML input file. ";
      std::cerr << boost::typeindex::type_id<T>() << " type was expected. \n";
      abort();
    }
    return;
  }

  template <class T> void read_verb(T &x, const YAML::Node &nd, const std::string &name) {
    read(x, nd, name);
    std::cout << "## " << name << "=" << x << "\n";
    return;
  }

  template <class T>
  void read_optional(T &x, const YAML::Node &nd, const std::string &name) {
    if (nd[name]) {
      read<T>(x, nd, name);
    }
    return;
  }

  template <class T>
  void read_opt_verb(T &x, const YAML::Node &nd, const std::string &name) {
    if (nd[name]) {
      read<T>(x, nd, name);
    } else {
      std::cout << "## " << name << "=" << x << " (default)\n";
    }
    return;
  }

} // namespace YAML_parsing

namespace input_file_parsing {

  int validate_geometry(gp::general &gparams) {
    if (gparams.ndims > 4 || gparams.ndims < 2) {
      std::cerr << "2 <= ndims <= 4!" << std::endl;
      return 1;
    }
    if (gparams.Lx < 1 || gparams.Ly < 1 || gparams.Lz < 1 || gparams.Lt < 1) {
      std::cerr << "All box extends must be > 1!" << std::endl;
      return 1;
    }
    if (gparams.ndims == 2) { // flattening 'y' and 'z' directions
      gparams.Ly = 1;
      gparams.Lz = 1;
    }
    if (gparams.ndims == 3) { // flattening 'z' direction
      gparams.Lz = 1;
    }

    return 0;
  }

  namespace u1 {

    namespace hmc {
      int parse_input_file(const std::string &file,
                           gp::general &gparams,
                           gp::hmc_u1 &hparams) {
        const YAML::Node nd = YAML::LoadFile(file);

        namespace Yp = YAML_parsing;
        std::cout << "## Parsing input file: " << file << "\n";

        // general parameters
        Yp::read_verb<size_t>(gparams.Lx, nd["geometry"], "X");
        Yp::read_verb<size_t>(gparams.Ly, nd["geometry"], "Y");
        Yp::read_verb<size_t>(gparams.Lz, nd["geometry"], "Z");
        Yp::read_verb<size_t>(gparams.Lt, nd["geometry"], "T");
        Yp::read_verb<size_t>(gparams.ndims, nd["geometry"], "ndims");

        int gerr = validate_geometry(gparams);
        if (gerr > 0) {
          return gerr;
        }

        Yp::read_verb<double>(gparams.beta, nd["action"], "beta");
        Yp::read_verb<double>(gparams.m0, nd["action"], "mass");

        Yp::read_verb<size_t>(hparams.N_save, nd["hmc"], "nsave");
        Yp::read_verb<size_t>(hparams.N_meas, nd["hmc"], "nmeas");
        Yp::read_verb<size_t>(hparams.icounter, nd["hmc"], "counter");
        Yp::read_verb<size_t>(hparams.seed, nd["hmc"], "seed");
        Yp::read_verb<double>(hparams.heat, nd["hmc"], "heat");
        Yp::read_verb<std::string>(hparams.configfilename, nd["hmc"], "configname");

        // hmc-u1 parameters
        Yp::read_opt_verb<size_t>(hparams.N_rev, nd["hmc"], "N_rev");
        Yp::read_opt_verb<size_t>(hparams.n_steps, nd["hmc"], "n_steps");
        Yp::read_opt_verb<double>(hparams.tau, nd["hmc"], "tau");
        Yp::read_opt_verb<size_t>(hparams.exponent, nd["hmc"], "exponent");
        Yp::read_opt_verb<size_t>(hparams.integs, nd["hmc"], "integs");
        Yp::read_opt_verb<bool>(hparams.no_fermions, nd["hmc"], "no_fermions");
        Yp::read_opt_verb<std::string>(hparams.solver, nd["hmc"], "solver");
        Yp::read_opt_verb<double>(hparams.tolerance_cg, nd["hmc"], "tolerance_cg");
        Yp::read_opt_verb<size_t>(hparams.solver_verbosity, nd["hmc"],
                                  "solver_verbosity");
        Yp::read_opt_verb<size_t>(hparams.seed_pf, nd["hmc"], "seed_pf");
        Yp::read_opt_verb<std::string>(hparams.outdir, nd["hmc"], "outdir");

        return 0;
      }

    } // namespace hmc

    namespace measure {

      int parse_input_file(const std::string &file,
                           gp::general &gparams,
                           gp::measure_u1 &mparams) {
        const YAML::Node nd = YAML::LoadFile(file);

        namespace Yp = YAML_parsing;
        std::cout << "## Parsing input file: " << file << "\n";

        // general parameters
        Yp::read_verb<size_t>(gparams.Lx, nd["geometry"], "X");
        Yp::read_verb<size_t>(gparams.Ly, nd["geometry"], "Y");
        Yp::read_verb<size_t>(gparams.Lz, nd["geometry"], "Z");
        Yp::read_verb<size_t>(gparams.Lt, nd["geometry"], "T");
        Yp::read_verb<size_t>(gparams.ndims, nd["geometry"], "ndims");

        int gerr = validate_geometry(gparams);
        if (gerr > 0) {
          return gerr;
        }
        Yp::read_verb<double>(gparams.beta, nd["action"], "beta");

        // measure-u1 parameters
        Yp::read_opt_verb<size_t>(mparams.nmeas, nd["measure"], "nmeas");
        Yp::read_opt_verb<size_t>(mparams.nstep, nd["measure"], "nstep");
        Yp::read_opt_verb<bool>(mparams.Wloop, nd["measure"], "Wloop");
        Yp::read_opt_verb<bool>(mparams.gradient, nd["measure"], "gradient");
        if (mparams.gradient) {
          Yp::read_opt_verb<double>(mparams.tmax, nd["measure"], "tmax");
        }
        Yp::read_opt_verb<std::string>(mparams.confdir, nd["measure"], "confdir");

        return 0;
      }

    } // namespace measure

  } // namespace u1
} // namespace input_file_parsing