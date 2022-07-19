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

  void parse_command_line(int ac, char *av[], std::string &input_file) {
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
      std::abort();
    }

    if (input_file == "NONE") {
      std::cout << "Error: invalid command line parameters.\nCheck using the "
                   "\"-h/--help\" option."
                << "\n";
      std::abort();
    }

    return;
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
    if (pparams.ndims < 4) {
      std::cerr << "## Warning: ndims==" << pparams.ndims << " --> flattening the \'z\' ";
      pparams.Lz = 1;
      std::string s_end = "";
      if (pparams.ndims < 3) {
        pparams.Ly = 1;
        std::cerr << "and \'y\' ";
        s_end = "s";
      }
      std::cerr << "direction" << s_end << "\n";
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
      std::abort();
    }
    return;
  }

  namespace u1 {
    namespace Yp = YAML_parsing;

    void parse_geometry(Yp::inspect_node &in, gp::physics &pparams) {
      YAML::Node R = in.get_root();

      const bool spec_L = (bool)R["geometry"]["L"];
      const bool spec_Lxyz =
        (R["geometry"]["X"] || R["geometry"]["Y"] || R["geometry"]["Z"]);

      if (spec_L ^ spec_Lxyz) {
        // either L^3*T or Lx*Ly*Lz*T have been specified (but not both)

        if (spec_L) {
          size_t L;
          in.read_verb<size_t>(L, {"geometry", "L"});
          pparams.Lx = L;
          pparams.Ly = L;
          pparams.Lz = L;
        } else {
          in.read_verb<size_t>(pparams.Lx, {"geometry", "X"});
          in.read_verb<size_t>(pparams.Ly, {"geometry", "Y"});
          in.read_verb<size_t>(pparams.Lz, {"geometry", "Z"});
        }

      } else {
        std::cerr << "Error: check your input file. ";
        std::cerr
          << "Either you specify L=Lx=Ly=Lz or each dimension separately, not both.\n";
        std::cerr << "Aborting.\n";
        std::abort();
      }

      in.read_verb<size_t>(pparams.Lt, {"geometry", "T"});
      in.read_verb<size_t>(pparams.ndims, {"geometry", "ndims"});

      int gerr = validate_geometry(pparams);
      if (gerr > 0) {
        std::cerr << "Error: invalid geometry parameters. Check L (or X,Y,Z), T, ndims "
                     "in your input file.";
        std::cerr << "Aborting.\n";
        std::abort();
      }

      if (R["geometry"]["rotating_frame"]) {
        pparams.rotating_frame = true;
        pparams.flat_metric = false; // metric is not flat
        in.read_verb<double>(pparams.Omega, {"geometry", "rotating_frame", "Omega"});
      }

      return;
    }

    /**
     * @brief parsing parameters from node of glueball correlator
     * @param in input file parser
     * @param mgparams measure_glueball_parameters structure
     */
    void parse_glueball_measure(Yp::inspect_node &in,
                                const std::vector<std::string> &inner_tree,
                                gp::measure_glueball_u1 &mgparams) {
      in.set_InnerTree(inner_tree); // entering the glueball node

      in.read_verb<bool>(mgparams.doAPEsmear, {"do_APE_smearing"});
      if (mgparams.doAPEsmear) {
        in.read_verb<size_t>(mgparams.nAPEsmear, {"APE_smearing", "n"});
        in.read_verb<double>(mgparams.alphaAPEsmear, {"APE_smearing", "alpha"});
      }
      in.read_opt_verb<bool>(mgparams.lengthy_file_name, {"lenghty_file_name"});
      in.read_opt_verb<bool>(mgparams.use_res_dir, {"res_dir"});

      in.read_opt_verb<bool>(mgparams.loops_GEVP, {"interpolators", "loops_GEVP"});
      in.read_opt_verb<size_t>(mgparams.max_length_loops, {"max_length_loops"});

      in.read_opt_verb<bool>(mgparams.U_ij, {"interpolators", "U_ij"});
      in.read_opt_verb<bool>(mgparams.U_munu, {"interpolators", "U_munu"});

      in.set_InnerTree({}); // reset to previous state
    }

    namespace hmc {
      void parse_input_file(const std::string &file,
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
            if (nd["monomials"]["gauge"]["anisotropic"]) {
              pparams.anisotropic = true;
              in.read_opt_verb<double>(pparams.xi,
                                       {"monomials", "gauge", "anisotropic", "xi"});
            }
          }
          if (nd["operators"]) {
            // note: initializing the D*Ddag monomial makes sense only with a valid Dirac
            // operator
            if (nd["operators"]["staggered"] && nd["monomials"]["staggered_det_DDdag"]) {
              pparams.include_staggered_fermions = true;

              in.read_verb<double>(pparams.m0, {"operators", "staggered", "mass"});

              in.read_opt_verb<std::string>(
                hparams.solver, {"monomials", "staggered_det_DDdag", "solver"});
              in.read_opt_verb<double>(
                hparams.tolerance_cg,
                {"monomials", "staggered_det_DDdag", "tolerance_cg"});
              in.read_opt_verb<size_t>(
                hparams.solver_verbosity,
                {"monomials", "staggered_det_DDdag", "solver_verbosity"});
              in.read_opt_verb<size_t>(hparams.seed_pf,
                                       {"monomials", "staggered_det_DDdag", "seed_pf"});
            }
          }
        }

        // hmc-u1 parameters
        in.read_verb<size_t>(hparams.N_save, {"hmc", "n_save"});
        in.read_verb<size_t>(hparams.n_meas, {"hmc", "n_meas"});

        in.read_opt_verb<bool>(hparams.do_hmc, {"hmc", "do_hmc"});

        if (nd["hmc"]["restart"] && nd["hmc"]["heat"]) {
          std::cerr << "Error: "
                    << "'restart' and 'heat' conditions are incompatible in the hmc. "
                    << "Aborting.\n";
          abort();
        }
        if (!nd["hmc"]["restart"] && !nd["hmc"]["heat"]) {
          std::cerr << "Error: "
                    << "Please pass either 'restart' or 'heat' to the hmc. "
                    << "Aborting.\n";
          abort();
        }

        in.read_opt_verb<bool>(hparams.restart, {"hmc", "restart"});
        in.read_opt_verb<bool>(hparams.heat, {"hmc", "heat"});

        in.read_opt_verb<size_t>(hparams.seed, {"hmc", "seed"});
        in.read_opt_verb<std::string>(hparams.configfilename, {"hmc", "configname"});
        in.read_opt_verb<std::string>(hparams.conf_dir, {"hmc", "conf_dir"});
        in.read_opt_verb<std::string>(hparams.conf_basename, {"hmc", "conf_basename"});
        in.read_opt_verb<bool>(hparams.lenghty_conf_name, {"hmc", "lenghty_conf_name"});

        in.read_opt_verb<size_t>(hparams.beta_str_width, {"hmc", "beta_str_width"});
        validate_beta_str_width(hparams.beta_str_width);

        // integrator parameters
        in.read_opt_verb<size_t>(hparams.N_rev, {"integrator", "N_rev"});
        in.read_opt_verb<size_t>(hparams.n_steps, {"integrator", "n_steps"});
        in.read_opt_verb<double>(hparams.tau, {"integrator", "tau"});
        in.read_opt_verb<size_t>(hparams.exponent, {"integrator", "exponent"});
        in.read_opt_verb<std::string>(hparams.integrator, {"integrator", "name"});

        if (nd["omeas"]) {
          hparams.make_omeas = true;

          hparams.omeas.res_dir = hparams.conf_dir; // default
          in.read_opt_verb<std::string>(hparams.omeas.res_dir, {"omeas", "res_dir"});

          in.read_opt_verb<size_t>(hparams.omeas.verbosity, {"omeas", "verbosity"});
          in.read_opt_verb<size_t>(hparams.omeas.icounter, {"omeas", "icounter"});
          in.read_opt_verb<size_t>(hparams.omeas.nstep, {"omeas", "nstep"});

          if (nd["omeas"]["pion_staggered"]) {
            hparams.omeas.pion_staggered = true;
            in.read_verb<double>(hparams.omeas.m0, {"omeas", "pion_staggered", "mass"});
          }

          in.read_opt_verb<bool>(hparams.omeas.Wloop, {"omeas", "Wloop"});

          if (nd["omeas"]["glueball"]) {
            hparams.omeas.glueball.do_measure = true;
            parse_glueball_measure(in, {"omeas", "glueball"}, hparams.omeas.glueball);
          }

          if (nd["omeas"]["gradient_flow"]) {
            hparams.omeas.gradient_flow = true;
            in.read_opt_verb<double>(hparams.omeas.epsilon_gradient_flow,
                                     {"omeas", "gradient_flow", "epsilon"});
            in.read_verb<double>(hparams.omeas.tmax, {"omeas", "gradient_flow", "tmax"});
          }
        }

        in.finalize();
        return;
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
          in.read_opt_verb<double>(pparams.xi,
                                   {"monomials", "gauge", "anisotropic", "xi"});
        }

        // measure-u1 parameters
        in.read_opt_verb<size_t>(mparams.n_meas, {"omeas", "n_meas"});
        in.read_opt_verb<size_t>(mparams.nstep, {"omeas", "nstep"});
        in.read_opt_verb<size_t>(mparams.icounter, {"omeas", "icounter"});
        in.read_opt_verb<size_t>(mparams.seed, {"omeas", "seed"});
        in.read_opt_verb<bool>(mparams.Wloop, {"omeas", "Wloop"});
        // optional parameters for gradient
        if (nd["omeas"]["gradient_flow"]) {
          mparams.gradient_flow = true;
          in.read_opt_verb<double>(mparams.epsilon_gradient_flow,
                                   {"omeas", "gradient_flow", "epsilon"});

          in.read_opt_verb<double>(mparams.tmax, {"omeas", "gradient_flow", "tmax"});
        }
        // optional parameters for pion
        if (nd["omeas"]["pion_staggered"]) {
          pparams.include_staggered_fermions = true;

          mparams.pion_staggered = true;
          in.read_verb<double>(pparams.m0, {"omeas", "pion_staggered", "mass"});

          in.read_opt_verb<std::string>(mparams.solver,
                                        {"omeas", "pion_staggered", "solver"});
          in.read_opt_verb<double>(mparams.tolerance_cg,
                                   {"omeas", "pion_staggered", "tolerance_cg"});
          in.read_opt_verb<size_t>(mparams.solver_verbosity,
                                   {"omeas", "pion_staggered", "solver_verbosity"});
          in.read_opt_verb<size_t>(mparams.seed_pf,
                                   {"omeas", "pion_staggered", "seed_pf"});
        }

        if (nd["omeas"]["glueball"]) {
          mparams.glueball.do_measure = true;
          parse_glueball_measure(in, {"omeas", "glueball"}, mparams.glueball);
        }

        // optional parameters for potentials
        if (nd["omeas"]["potential"]) {
          in.read_opt_verb<bool>(mparams.potential, {"omeas", "potential", "potential"});
          in.read_opt_verb<bool>(mparams.potentialsmall,
                                 {"omeas", "potential", "potentialsmall"});
          in.read_opt_verb<bool>(mparams.append, {"omeas", "potential", "append"});
          in.read_opt_verb<bool>(mparams.smear_spatial_only,
                                 {"omeas", "potential", "smear_spatial_only"});
          in.read_opt_verb<bool>(mparams.smear_temporal_only,
                                 {"omeas", "potential", "smear_temporal_only"});
          in.read_opt_verb<size_t>(mparams.n_apesmear,
                                   {"omeas", "potential", "n_apesmear"});
          in.read_opt_verb<double>(mparams.alpha, {"omeas", "potential", "alpha"});
          in.read_opt_verb<double>(mparams.sizeWloops,
                                   {"omeas", "potential", "sizeWloops"});
        }

        in.read_opt_verb<std::string>(mparams.conf_dir, {"omeas", "conf_dir"});
        in.read_opt_verb<std::string>(mparams.res_dir, {"omeas", "res_dir"});

        in.read_opt_verb<std::string>(mparams.conf_basename, {"omeas", "conf_basename"});

        in.read_opt_verb<size_t>(mparams.beta_str_width, {"omeas", "beta_str_width"});
        validate_beta_str_width(mparams.beta_str_width);

        in.finalize();

        in.finalize();
        return 0;
      }

    } // namespace measure

    namespace metropolis {
      /**
       * @brief check if N_hit >= 1, otherwise it aborts
       * @param n N_hit (number of times an update per link is attempted)
       */
      void validate_N_hit(const size_t &n) {
        if (n < 1) {
          std::cerr << "Error: N_hit should be at least 1, otherwise nothing will happen";
          std::cerr << "Aborting.\n";
          std::abort();
        }
        return;
      }

      int parse_input_file(const std::string &file,
                           gp::physics &pparams,
                           gp::metropolis_u1 &mcparams) {
        std::cout << "## Parsing input file: " << file << "\n";
        const YAML::Node nd = YAML::LoadFile(file);
        Yp::inspect_node in(nd);

        parse_geometry(in, pparams);

        // beta, xi value from the gauge action
        in.read_verb<double>(pparams.beta, {"monomials", "gauge", "beta"});
        in.read_opt_verb<bool>(pparams.anisotropic,
                               {"monomials", "gauge", "anisotropic"});
        if (pparams.anisotropic) {
          in.read_opt_verb<double>(pparams.xi, {"monomials", "gauge", "xi"});
        }

        // measure-u1 parameters
        in.read_opt_verb<size_t>(mcparams.n_meas, {"metropolis", "n_meas"});
        in.read_opt_verb<size_t>(mcparams.N_save, {"metropolis", "N_save"});
        in.read_opt_verb<size_t>(mcparams.seed, {"metropolis", "seed"});

        in.read_opt_verb<std::string>(mcparams.conf_dir, {"metropolis", "conf_dir"});
        in.read_opt_verb<std::string>(mcparams.conf_basename,
                                      {"metropolis", "conf_basename"});
        in.read_opt_verb<bool>(mcparams.lenghty_conf_name,
                               {"metropolis", "lenghty_conf_name"});
        in.read_opt_verb<size_t>(mcparams.beta_str_width,
                                 {"metropolis", "beta_str_width"});
        validate_beta_str_width(mcparams.beta_str_width);
        in.read_opt_verb<bool>(mcparams.restart, {"metropolis", "restart"});
        if (mcparams.restart) {
          in.read_verb<size_t>(mcparams.icounter, {"metropolis", "icounter"});
        } else {
          in.read_verb<double>(mcparams.heat, {"metropolis", "heat"});
        }

        in.read_verb<double>(mcparams.delta, {"metropolis", "delta"});
        in.read_opt_verb<size_t>(mcparams.N_hit, {"metropolis", "N_hit"});
        validate_N_hit(mcparams.N_hit);

        in.finalize();

        return 0;
      }

    } // namespace metropolis

  } // namespace u1
} // namespace input_file_parsing
