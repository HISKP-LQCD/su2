// parse_input_file.cc

#include <iostream>
#include <string>

#include "parse_input_file.hh"

namespace input_file_parsing {

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

        if (nd["begin_monomials"]) {
          if (nd["begin_monomials"]["gauge"]) {
            pparams.include_gauge = true;
            in.read_verb<double>(pparams.beta, {"begin_monomials", "gauge", "beta"});
          }
          // note: initializing the D*Ddag monomial makes sense only with a valid Dirac
          // operator
          if (nd["begin_operators"]["staggered"] && nd["begin_monomials"]["staggered_det_DDdag"]) {
            pparams.include_staggered_fermions = true;

            in.read_verb<double>(pparams.m0, {"begin_operators", "staggered", "mass"});

            in.read_opt_verb<std::string>(hparams.solver, {"begin_monomials", "staggered_det_DDdag",
                                           "solver"});
            in.read_opt_verb<double>(hparams.tolerance_cg, {"begin_monomials", "staggered_det_DDdag",
                                      "tolerance_cg"});
            in.read_opt_verb<size_t>(hparams.solver_verbosity,
                                      {"begin_monomials", "staggered_det_DDdag", "solver_verbosity"});
            in.read_opt_verb<size_t>(hparams.seed_pf, {"begin_monomials", "staggered_det_DDdag",
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
        in.read_opt_verb<size_t>(hparams.N_rev, {"begin_integrator", "N_rev"});
        in.read_opt_verb<size_t>(hparams.n_steps, {"begin_integrator", "n_steps"});
        in.read_opt_verb<double>(hparams.tau, {"begin_integrator", "tau"});
        in.read_opt_verb<size_t>(hparams.exponent, {"begin_integrator", "exponent"});
        in.read_opt_verb<std::string>(hparams.integrator, {"begin_integrator",
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
        in.read_verb<double>(pparams.beta, {"begin_monomials", "gauge", "beta"});
        in.read_opt_verb<bool>(pparams.anisotropic, {"begin_monomials", "gauge",
                                "anisotropic"});
        if (pparams.anisotropic) {
          in.read_opt_verb<double>(pparams.xi, {"begin_monomials", "gauge", "xi"});
        }

        // measure-u1 parameters
        in.read_opt_verb<size_t>(mparams.n_meas, {"begin_measurements", "n_meas"});
        in.read_opt_verb<size_t>(mparams.nstep, {"begin_measurements", "nstep"});
        in.read_opt_verb<size_t>(mparams.icounter, {"begin_measurements", "icounter"});
        in.read_opt_verb<size_t>(mparams.seed, {"begin_measurements", "seed"});
        in.read_opt_verb<bool>(mparams.Wloop, {"begin_measurements", "Wloop"});
        // optional parameters for gradient
        if (nd["begin_measurements"]["gradient"]) {
          in.read_opt_verb<double>(mparams.tmax, {"begin_measurements", "gradient", "tmax"});
        }
        // optional parameters for potentials
        if (nd["begin_measurements"]["potential"]) {
          in.read_opt_verb<bool>(mparams.potential, {"begin_measurements", "potential", "potential"});
          in.read_opt_verb<bool>(mparams.potentialsmall, {"begin_measurements", "potential", "potentialsmall"});
          in.read_opt_verb<bool>(mparams.append, {"begin_measurements", "potential", "append"});
          in.read_opt_verb<bool>(mparams.smear_spatial_only, {"begin_measurements", "potential",
                                  "smear_spatial_only"});
          in.read_opt_verb<size_t>(mparams.n_apesmear, {"begin_measurements", "potential", "n_apesmear"});
          in.read_opt_verb<double>(mparams.alpha, {"begin_measurements", "potential", "alpha"});
          in.read_opt_verb<double>(mparams.sizeWloops, {"begin_measurements", "potential", "sizeWloops"});
        }

        in.read_opt_verb<std::string>(mparams.confdir, {"begin_measurements", "confdir"});
        in.read_opt_verb<std::string>(mparams.resdir, {"begin_measurements", "resdir"});

        in.read_opt_verb<std::string>(mparams.conf_basename, {"begin_measurements", "conf_basename"});

        in.read_opt_verb<size_t>(mparams.beta_str_width, {"begin_measurements", "beta_str_width"});
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
        in.read_verb<double>(pparams.beta, {"begin_monomials", "gauge", "beta"});
        in.read_opt_verb<bool>(pparams.anisotropic, {"begin_monomials", "gauge","anisotropic"});
        if (pparams.anisotropic) {
          in.read_opt_verb<double>(pparams.xi, {"begin_monomials", "gauge", "xi"});
        }

        // measure-u1 parameters
        in.read_opt_verb<size_t>(mcparams.n_meas, {"begin_metropolis", "n_meas"});
        in.read_opt_verb<size_t>(mcparams.N_save, {"begin_metropolis", "N_save"});
        in.read_opt_verb<size_t>(mcparams.icounter, {"begin_metropolis", "icounter"});
        in.read_opt_verb<size_t>(mcparams.seed, {"begin_metropolis", "seed"});

        in.read_opt_verb<std::string>(mcparams.outdir, {"begin_metropolis", "outdir"});
        in.read_opt_verb<std::string>(mcparams.conf_basename, {"begin_metropolis", "conf_basename"});
        in.read_opt_verb<size_t>(mcparams.beta_str_width, {"begin_metropolis", "beta_str_width"});
        in.read_opt_verb<bool>(mcparams.restart, {"begin_metropolis", "restart"});
        if (mcparams.restart) {
          in.read_opt_verb<std::string>(mcparams.configfilename, {"begin_metropolis", "configfilename"});
        }

        in.read_verb<double>(mcparams.heat, {"begin_metropolis", "heat"});
        in.read_verb<double>(mcparams.delta, {"begin_metropolis", "delta"});
        in.read_opt_verb<size_t>(mcparams.N_hit, {"begin_metropolis", "N_hit"});
        
        in.finalize();

        return 0;
      }

    } // namespace metropolis

  } // namespace u1
} // namespace input_file_parsing