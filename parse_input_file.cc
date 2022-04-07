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

  int validate_geometry(gp::physics &pparams) {
    if (pparams.ndims > 4 || pparams.ndims < 2) {
      std::cerr << "2 <= ndims <= 4!" << std::endl;
      return 1;
    }
    if (pparams.Lx < 1 || pparams.Ly < 1 || pparams.Lz < 1 || pparams.Lt < 1) {
      std::cerr << "All box extends must be > 1!" << std::endl;
      return 1;
    }
    if (pparams.ndims == 2) { // flattening 'y' and 'z' directions
      pparams.Ly = 1;
      pparams.Lz = 1;
    }
    if (pparams.ndims == 3) { // flattening 'z' direction
      pparams.Lz = 1;
    }

    return 0;
  }

  namespace u1 {
    namespace Yp = YAML_parsing;

    int parse_geometry(const YAML::Node &nd, gp::physics &pparams){

        // physics parameters
        Yp::read_verb<size_t>(pparams.Lx, nd["geometry"], "X");
        Yp::read_verb<size_t>(pparams.Ly, nd["geometry"], "Y");
        Yp::read_verb<size_t>(pparams.Lz, nd["geometry"], "Z");
        Yp::read_verb<size_t>(pparams.Lt, nd["geometry"], "T");
        Yp::read_verb<size_t>(pparams.ndims, nd["geometry"], "ndims");

        int gerr = validate_geometry(pparams);
        if (gerr > 0) {
          return gerr;
        }

        return gerr;
    }

    namespace hmc {
      int parse_input_file(const std::string &file,
                           gp::physics &pparams,
                           gp::hmc_u1 &hparams) {
        std::cout << "## Parsing input file: " << file << "\n";
        const YAML::Node nd = YAML::LoadFile(file);

        parse_geometry(nd, pparams);

        if (nd["begin_monomials"]) {
          const YAML::Node &nBM = nd["begin_monomials"];
          if (nBM["gauge"]) {
            pparams.include_gauge = true;
            Yp::read_verb<double>(pparams.beta, nBM["gauge"], "beta");
          }
          // note: initializing the D*Ddag monomial makes sense only with a valid Dirac
          // operator
          if (nd["begin_operators"]["staggered"] && nBM["staggered_det_DDdag"]) {
            pparams.include_staggered_fermions = true;
            const YAML::Node &nBO = nd["begin_operators"];

            Yp::read_verb<double>(pparams.m0, nBO["staggered"], "mass");

            Yp::read_opt_verb<std::string>(hparams.solver, nBM["staggered_det_DDdag"],
                                           "solver");
            Yp::read_opt_verb<double>(hparams.tolerance_cg, nBM["staggered_det_DDdag"],
                                      "tolerance_cg");
            Yp::read_opt_verb<size_t>(hparams.solver_verbosity,
                                      nBM["staggered_det_DDdag"], "solver_verbosity");
            Yp::read_opt_verb<size_t>(hparams.seed_pf, nBM["staggered_det_DDdag"],
                                      "seed_pf");
          }
        }

        // hmc-u1 parameters
        Yp::read_verb<size_t>(hparams.N_save, nd["hmc"], "n_save");
        Yp::read_verb<size_t>(hparams.N_meas, nd["hmc"], "n_meas");
        Yp::read_verb<size_t>(hparams.icounter, nd["hmc"], "counter");
        Yp::read_verb<double>(hparams.heat, nd["hmc"], "heat");
        Yp::read_opt_verb<size_t>(hparams.seed, nd["hmc"], "seed");
        Yp::read_opt_verb<std::string>(hparams.configfilename, nd["hmc"], "configname");
        Yp::read_opt_verb<std::string>(hparams.outdir, nd["hmc"], "outdir");
        Yp::read_opt_verb<std::string>(hparams.conf_basename, nd["hmc"], "conf_basename");
        Yp::read_opt_verb<size_t>(hparams.beta_str_width, nd["hmc"], "beta_str_width");

        // integrator parameters
        Yp::read_opt_verb<size_t>(hparams.N_rev, nd["begin_integrator"], "N_rev");
        Yp::read_opt_verb<size_t>(hparams.n_steps, nd["begin_integrator"], "n_steps");
        Yp::read_opt_verb<double>(hparams.tau, nd["begin_integrator"], "tau");
        Yp::read_opt_verb<size_t>(hparams.exponent, nd["begin_integrator"], "exponent");
        Yp::read_opt_verb<std::string>(hparams.integrator, nd["begin_integrator"],
                                       "name");


        return 0;
      }

    } // namespace hmc

    namespace measure {

      int parse_input_file(const std::string &file,
                           gp::physics &pparams,
                           gp::measure_u1 &mparams) {
        std::cout << "## Parsing input file: " << file << "\n";
        const YAML::Node nd = YAML::LoadFile(file);

        parse_geometry(nd, pparams);


        // beta, xi value from the gauge action
        Yp::read_verb<double>(pparams.beta, nd["begin_monomials"]["gauge"], "beta");
        Yp::read_verb<double>(pparams.xi, nd["begin_monomials"]["gauge"], "xi");
        Yp::read_verb<bool>(pparams.anisotropic, nd["begin_monomials"]["gauge"], "anisotropic");

        // measure-u1 parameters
        const YAML::Node &nMS = nd["begin_measurements"];
        Yp::read_opt_verb<size_t>(mparams.nmeas, nMS, "nmeas");
        Yp::read_opt_verb<size_t>(mparams.nstep, nMS, "nstep");
        Yp::read_opt_verb<size_t>(mparams.icounter, nMS, "icounter");
        Yp::read_opt_verb<bool>(mparams.Wloop, nMS, "Wloop");
        //optional parameters for gradient
        if (nMS["gradient"]) {
          Yp::read_opt_verb<double>(mparams.tmax, nMS["gradient"], "tmax");
        }
        //optional parameters for potentials
        if (nMS["potential"]){
          Yp::read_opt_verb<bool>(mparams.potential, nMS["potential"], "potential");
          Yp::read_opt_verb<bool>(mparams.potentialsmall, nMS["potential"], "potentialsmall");
          Yp::read_opt_verb<bool>(mparams.append, nMS["potential"], "append");
          Yp::read_opt_verb<bool>(mparams.smearspacial, nMS["potential"], "smearspacial");
          Yp::read_opt_verb<size_t>(mparams.n_apesmear, nMS["potential"], "n_apesmear");
          Yp::read_opt_verb<double>(mparams.alpha, nMS["potential"], "alpha");
          Yp::read_opt_verb<double>(mparams.sizeloops, nMS["potential"], "sizeloops");
        }
        
        Yp::read_opt_verb<std::string>(mparams.confdir, nMS, "confdir");

        Yp::read_opt_verb<std::string>(mparams.conf_basename, nMS, "conf_basename");
        Yp::read_opt_verb<size_t>(mparams.beta_str_width, nMS, "beta_str_width");

        return 0;
      }

    } // namespace measure

  } // namespace u1
} // namespace input_file_parsing
