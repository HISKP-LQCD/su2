// parse_input_file.cc
/**
 * @file parse_input_file.cc
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief definitions of parse_input_file.hh 
 * @version 0.1
 * @date 2022-05-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include <string>

#include "parse_input_file.hh"

namespace input_file_parsing {

  int parse_command_line(int ac, char *av[], std::string &input_file) {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "produce this help message")(
      "file,f", po::value<std::string>(&input_file)->default_value("NONE"),
      "yaml input file");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
    }

    if(input_file == "NONE"){
      std::cout << "Error: invalid command line parameters.\nCheck using the \"-h/--help\" option." << "\n";
      return 1;
    }

    return 0;
  }

  int validate_geometry(gp::physics &pparams) {
    if (pparams.ndims > 4 || pparams.ndims < 2) {
      std::cerr << "2 <= ndims <= 4!" << std::endl;
      return 1;
    }
    if (pparams.Lx < 1 || pparams.Ly < 1 || pparams.Lz < 1 || pparams.Lt < 1) {
      std::cerr << "All box extents must be > 1!" << std::endl;
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

  /**
   * @brief check if beta_str_width >= 1, otherwise it aborts
   * @param n width (number of digits after the dot)
   */
  void validate_beta_str_width(const size_t &n) {
    if (n < 1) {
      std::cerr << "Error: beta string width in the output file should have at least "
                   "decimal 1 digit.";
      std::cerr << "Aborting.\n";
      abort();
    }
    return;
  }

  namespace u1 {
    namespace Yp = YAML_parsing;

    void parse_geometry(Yp::inspect_node& in, gp::physics &pparams) {

      in.read_verb<size_t>(pparams.Lx, {"geometry", "X"});
      in.read_verb<size_t>(pparams.Ly, {"geometry", "Y"});
      in.read_verb<size_t>(pparams.Lz, {"geometry", "Z"});
      in.read_verb<size_t>(pparams.Lt, {"geometry", "T"});
      in.read_verb<size_t>(pparams.ndims, {"geometry", "ndims"});

      int gerr = validate_geometry(pparams);
      if (gerr > 0) {
        std::cerr
          << "Error: invalid geometry parameters. Check X,Y,Z,ndims in your input file.";
        std::cerr << "Aborting.\n";
        abort();
      }

      return;
    }

    namespace hmc {
      int parse_input_file(const std::string &file,
                           gp::physics &pparams,
                           gp::hmc_u1 &hparams) {
        std::cout << "## Parsing input file: " << file << "\n";
        const YAML::Node nd = YAML::LoadFile(file);
        Yp::inspect_node in(nd);

        parse_geometry(in, pparams);

        if (nd["monomials"]) {
          if (nd["monomials"]["gauge"]) {
            pparams.include_gauge = true;
            in.read_verb<double>(pparams.beta, {"monomials", "gauge", "beta"});
          }
          // note: initializing the D*Ddag monomial makes sense only with a valid Dirac
          // operator
          if (nd["operators"]["staggered"] && nd["monomials"]["staggered_det_DDdag"]) {
            pparams.include_staggered_fermions = true;

            in.read_verb<double>(pparams.m0, {"operators", "staggered", "mass"});

            in.read_opt_verb<std::string>(hparams.solver, {"monomials", "staggered_det_DDdag",
                                           "solver"});
            in.read_opt_verb<double>(hparams.tolerance_cg, {"monomials", "staggered_det_DDdag",
                                      "tolerance_cg"});
            in.read_opt_verb<size_t>(hparams.solver_verbosity,
                                      {"monomials", "staggered_det_DDdag", "solver_verbosity"});
            in.read_opt_verb<size_t>(hparams.seed_pf, {"monomials", "staggered_det_DDdag",
                                      "seed_pf"});
          }
        }

        // hmc-u1 parameters
        in.read_verb<size_t>(hparams.N_save, {"hmc", "n_save"});
        in.read_verb<size_t>(hparams.n_meas, {"hmc", "n_meas"});
        in.read_verb<size_t>(hparams.icounter, {"hmc", "counter"});
        in.read_verb<bool>(hparams.heat, {"hmc", "heat"});
        in.read_opt_verb<size_t>(hparams.seed, {"hmc", "seed"});
        in.read_opt_verb<std::string>(hparams.configfilename, {"hmc", "configname"});
        in.read_opt_verb<std::string>(hparams.outdir, {"hmc", "outdir"});
        in.read_opt_verb<std::string>(hparams.conf_basename, {"hmc", "conf_basename"});

        in.read_opt_verb<size_t>(hparams.beta_str_width, {"hmc", "beta_str_width"});
        validate_beta_str_width(hparams.beta_str_width);

        // integrator parameters
        in.read_opt_verb<size_t>(hparams.N_rev, {"integrator", "N_rev"});
        in.read_opt_verb<size_t>(hparams.n_steps, {"integrator", "n_steps"});
        in.read_opt_verb<double>(hparams.tau, {"integrator", "tau"});
        in.read_opt_verb<size_t>(hparams.exponent, {"integrator", "exponent"});
        in.read_opt_verb<std::string>(hparams.integrator, {"integrator",
                                       "name"});

        in.finalize();
        return 0;
      }

    } // namespace hmc

    namespace measure {

      int parse_input_file(const std::string &file,
                           gp::physics &pparams,
                           gp::measure_u1 &mparams) {
        std::cout << "## Parsing input file: " << file << "\n";
        const YAML::Node nd = YAML::LoadFile(file);
        Yp::inspect_node in(nd);

        parse_geometry(in, pparams);

        // beta, xi value from the gauge action
        in.read_verb<double>(pparams.beta, {"monomials", "gauge", "beta"});
        if (nd["monomials"]["gauge"]["anisotropic"]) {
          pparams.anisotropic = true;
          in.read_opt_verb<double>(pparams.xi, {"monomials", "gauge", "anisotropic", "xi"});
        }

        // measure-u1 parameters
        in.read_opt_verb<size_t>(mparams.n_meas, {"measurements", "n_meas"});
        in.read_opt_verb<size_t>(mparams.nstep, {"measurements", "nstep"});
        in.read_opt_verb<size_t>(mparams.icounter, {"measurements", "icounter"});
        in.read_opt_verb<size_t>(mparams.seed, {"measurements", "seed"});
        in.read_opt_verb<bool>(mparams.Wloop, {"measurements", "Wloop"});
        // optional parameters for gradient
        if (nd["measurements"]["gradient"]) {
          in.read_opt_verb<double>(mparams.tmax, {"measurements", "gradient", "tmax"});
        }
        // optional parameters for potentials
        if (nd["measurements"]["potential"]) {
          mparams.potential=true;
          in.read_opt_verb<bool>(mparams.potentialsmall, {"measurements", "potential", "potentialsmall"});
          in.read_opt_verb<bool>(mparams.append, {"measurements", "potential", "append"});
          in.read_opt_verb<bool>(mparams.smear_spatial_only, {"measurements", "potential",
                                  "smear_spatial_only"});
          in.read_opt_verb<size_t>(mparams.n_apesmear, {"measurements", "potential", "n_apesmear"});
          in.read_opt_verb<double>(mparams.alpha, {"measurements", "potential", "alpha"});
          in.read_opt_verb<double>(mparams.sizeWloops, {"measurements", "potential", "sizeWloops"});
        }

        in.read_opt_verb<std::string>(mparams.confdir, {"measurements", "confdir"});
        in.read_opt_verb<std::string>(mparams.resdir, {"measurements", "resdir"});

        in.read_opt_verb<std::string>(mparams.conf_basename, {"measurements", "conf_basename"});

        in.read_opt_verb<size_t>(mparams.beta_str_width, {"measurements", "beta_str_width"});
        validate_beta_str_width(mparams.beta_str_width);

        in.finalize();
        return 0;
      }

    } // namespace measure

    namespace metropolis {

      int parse_input_file(const std::string &file,
                           gp::physics &pparams,
                           gp::metropolis_u1 &mcparams) {
        std::cout << "## Parsing input file: " << file << "\n";
        const YAML::Node nd = YAML::LoadFile(file);
        Yp::inspect_node in(nd);

        parse_geometry(in, pparams);

        // beta, xi value from the gauge action
        in.read_verb<double>(pparams.beta, {"monomials", "gauge", "beta"});
        in.read_opt_verb<bool>(pparams.anisotropic, {"monomials", "gauge","anisotropic"});
        if (pparams.anisotropic) {
          in.read_opt_verb<double>(pparams.xi, {"monomials", "gauge", "xi"});
        }

        // measure-u1 parameters
        in.read_opt_verb<size_t>(mcparams.n_meas, {"metropolis", "n_meas"});
        in.read_opt_verb<size_t>(mcparams.N_save, {"metropolis", "N_save"});
        in.read_opt_verb<size_t>(mcparams.icounter, {"metropolis", "icounter"});
        in.read_opt_verb<size_t>(mcparams.seed, {"metropolis", "seed"});

        in.read_opt_verb<std::string>(mcparams.outdir, {"metropolis", "outdir"});
        in.read_opt_verb<std::string>(mcparams.conf_basename, {"metropolis", "conf_basename"});
        in.read_opt_verb<size_t>(mcparams.beta_str_width, {"metropolis", "beta_str_width"});
        in.read_opt_verb<bool>(mcparams.restart, {"metropolis", "restart"});
        if (mcparams.restart) {
          in.read_opt_verb<std::string>(mcparams.configfilename, {"metropolis", "configfilename"});
        }

        in.read_verb<double>(mcparams.heat, {"metropolis", "heat"});
        in.read_verb<double>(mcparams.delta, {"metropolis", "delta"});
        in.read_opt_verb<size_t>(mcparams.N_hit, {"metropolis", "N_hit"});
        
        in.finalize();

        return 0;
      }

    } // namespace metropolis

  } // namespace u1
} // namespace input_file_parsing