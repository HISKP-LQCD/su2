/**
 * @file hmc-u1.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief class for the hmc algorithm for the U(1) gauge theory
 * @version 0.1
 * @date 2022-09-02
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "md_update.hh"
#include "program-u1.hpp"

#include "rotating-gaugemonomial.hpp"

#include "detDDdag_monomial.hh"

namespace u1 {

  class hmc_algo : public program<gp::hmc_u1> {
  private:
    std::string conf_path_basename;

    bool g_heat; // hot or cold starting configuration
    size_t g_icounter = 0; // 1st configuration(trajectory) to load from
    double normalisation;

    double rate = 0.; // acceptance rate

    std::list<monomial<double, _u1> *> monomial_list;

    flat_spacetime::gaugemonomial<double, _u1> *gm;
    rotating_spacetime::gauge_monomial<double, _u1> *gm_rot;

    kineticmonomial<double, _u1> *km;

    staggered::detDDdag_monomial<double, _u1> *detDDdag;

    integrator<double, _u1> *md_integ;
    md_params mdparams; // Molecular Dynamics parameters

  public:
    hmc_algo() {}
    ~hmc_algo() {
      free(gm);
      free(gm_rot);
      free(km);
      free(detDDdag);
      free(md_integ);
    }

    void print_program_info() const {
      std::cout << "## HMC Algorithm for U(1) gauge theory\n";
    }

    void parse_input_file() {
      namespace in_hmc = input_file_parsing::u1::hmc;
      in_hmc::parse_input_file(input_file, pparams, sparams);
      conf_path_basename = io::get_conf_path_basename(pparams, sparams);
      (*this).omeas = (*this).sparams.omeas;
    }

    /**
     * @brief initialize the potential filenames attributes
     */
    void set_potential_filenames() {
      // filename needed for saving results from potential and potentialsmall
      filename_fine = io::measure::get_filename_fine(pparams, sparams.omeas);
      filename_coarse = io::measure::get_filename_coarse(pparams, sparams.omeas);
      filename_nonplanar = io::measure::get_filename_nonplanar(pparams, sparams.omeas);

      // write explanatory headers into result-files, also check if measuring routine is
      // implemented for given dimension
      if (sparams.omeas.potentialplanar) {
        io::measure::set_header_planar(pparams, sparams.omeas, filename_coarse,
                                       filename_fine);
      }
      if (sparams.omeas.potentialnonplanar) {
        io::measure::set_header_nonplanar(pparams, sparams.omeas, filename_nonplanar);
      }
      return;
    }

    template <class Group, class URNG>
    std::vector<double> sweep(const gp::physics &pparams,
                              gaugeconfig<Group> &U,
                              std::vector<URNG> engines,
                              const double &delta,
                              const size_t &N_hit,
                              const double &beta,
                              const double &xi = 1.0,
                              const bool &anisotropic = false) {
      if (pparams.flat_metric) {
        return flat_spacetime::sweep(U, engines, delta, N_hit, pparams.beta, pparams.xi,
                                     pparams.anisotropic);
      }
      if (pparams.rotating_frame) {
        return rotating_spacetime::sweep(U, pparams.Omega, engines, delta, N_hit,
                                         pparams.beta, pparams.xi, pparams.anisotropic);
      } else {
        spacetime_lattice::fatal_error("Invalid metric when calling: ", __func__);
        return {};
      }
    }

    void init_gauge_conf() {
      gaugeconfig<_u1> U0(pparams.Lx, pparams.Ly, pparams.Lz, pparams.Lt, pparams.ndims,
                          pparams.beta);
      U = U0;

      if (sparams.restart) {
        std::cout << "## restart " << sparams.restart << "\n";
        const std::vector<std::string> v_ncc =
          io::hmc::read_nconf_counter(sparams.conf_dir);
        g_heat = boost::lexical_cast<bool>(v_ncc[0]);
        g_icounter = std::stoi(v_ncc[1]);
        std::string config_path = v_ncc[2];

        const size_t err = U.load(config_path);
        if (err != 0) {
          std::cout
            << "Error: failed to load initial gauge configuration for hmc. Aborting.\n";
          std::abort();
        }
      } else {
        g_heat = (sparams.heat == true) ? 1.0 : 0.0;
        g_icounter = 0;
        hotstart(U, sparams.seed, g_heat);
      }

      double plaquette = flat_spacetime::gauge_energy(U);
      double fac = 2. / U.getndims() / (U.getndims() - 1);
      normalisation = fac / U.getVolume() / double(U.getNc());

      std::cout << "## Normalization factor: A = 2/(d*(d-1)*N_lat*N_c) = "
                << std::scientific << std::setw(18) << std::setprecision(15)
                << normalisation << "\n";
      std::cout << "## Acceptance rate parcentage: rho = rate/(i+1)\n";

      std::cout << "## Initial Plaquette: " << plaquette * normalisation << std::endl;

      random_gauge_trafo(U, 654321);
      plaquette = flat_spacetime::gauge_energy(U);
      std::cout << "## Plaquette after rnd trafo: " << plaquette * normalisation
                << std::endl;
    }

    void open_output_data() {
      if (g_icounter == 0) {
        os.open(sparams.conf_dir + "/output.u1-hmc.data", std::ios::out);
      } else {
        os.open(sparams.conf_dir + "/output.u1-hmc.data", std::ios::app);
      }
    }

    // generate list of monomials
    void init_monomials() {
      gm = new flat_spacetime::gaugemonomial<double, _u1>(0, pparams.xi);
      gm_rot = new rotating_spacetime::gauge_monomial<double, _u1>(0, pparams.Omega);

      km = new kineticmonomial<double, _u1>(0);
      km->setmdpassive();
      monomial_list.push_back(km);

      detDDdag = new staggered::detDDdag_monomial<double, _u1>(
        0, pparams.m0, sparams.solver, sparams.tolerance_cg, sparams.seed_pf,
        sparams.solver_verbosity);

      if (pparams.include_gauge) {
        if (pparams.rotating_frame) {
          monomial_list.push_back(gm_rot);
        } else {
          monomial_list.push_back(gm);
        }
      }

      if (pparams.include_staggered_fermions) { // including S_F (fermionic) in the action
        monomial_list.push_back(detDDdag);
      }
    }

    void do_hmc_step(const int &i) {
      if (sparams.do_mcmc) {
        mdparams.disablerevtest();

        if (i > 0 && sparams.N_rev != 0 && (i) % sparams.N_rev == 0) {
          mdparams.enablerevtest();
        }
        // PRNG engine
        std::mt19937 engine(sparams.seed + i);
        // perform the MD update

        md_update(U, engine, mdparams, monomial_list, *md_integ);

        const double energy = flat_spacetime::gauge_energy(U);
        double E = 0., Q = 0.;
        flat_spacetime::energy_density(U, E, Q);
        rate += mdparams.getaccept();

        std::cout << i << " " << mdparams.getaccept() << " " << std::scientific
                  << std::setw(18) << std::setprecision(15) << energy * normalisation
                  << " " << std::setw(15) << mdparams.getdeltaH() << " " << std::setw(15)
                  << rate / static_cast<double>(i + 1) << " ";

        if (mdparams.getrevtest()) {
          std::cout << mdparams.getdeltadeltaH();
        } else {
          std::cout << "NA";
        }
        std::cout << " " << Q << std::endl;

        os << i << " " << mdparams.getaccept() << " " << std::scientific << std::setw(18)
           << std::setprecision(15) << energy * normalisation << " " << std::setw(15)
           << mdparams.getdeltaH() << " " << std::setw(15)
           << rate / static_cast<double>(i + 1) << " ";
        if (mdparams.getrevtest()) {
          os << mdparams.getdeltadeltaH();
        } else {
          os << "NA";
        }
        os << " " << Q << std::endl;
      }
    }

    void do_omeas_i(const size_t &i, const bool &do_omeas) {
      if (do_omeas) {
        if (sparams.omeas.Wloop) {
          if (sparams.omeas.verbosity > 0) {
            std::cout << "## online measuring: Wilson loop\n";
          }
          omeasurements::meas_wilson_loop<_u1>(U, i, sparams.conf_dir);
        }
        if (sparams.omeas.gradient_flow) {
          if (sparams.omeas.verbosity > 0) {
            std::cout << "## online measuring: Gradient flow\n";
          }
          omeasurements::meas_gradient_flow<_u1>(U, i, pparams, sparams.omeas);
        }

        if (sparams.omeas.pion_staggered) {
          if (sparams.omeas.verbosity > 0) {
            std::cout << "## online measuring: Pion correlator\n";
          }
          omeasurements::meas_pion_correlator<_u1>(U, i, pparams.m0, sparams.omeas);
        }

        if (sparams.omeas.glueball.do_measure) {
          if (sparams.omeas.verbosity > 0) {
            std::cout << "## online measuring: J^{PC} glueball correlators.\n";
          }
          if (sparams.omeas.glueball.loops_GEVP) {
            omeasurements::meas_glueball_correlator_GEVP<_u1>(U, i, sparams.omeas);
          }
          if (sparams.omeas.glueball.U_munu) {
            omeasurements::meas_glueball_correlator_U_munu<_u1>(U, i, sparams.omeas,
                                                                false);
          }
          if (sparams.omeas.glueball.U_ij) {
            omeasurements::meas_glueball_correlator_U_munu<_u1>(U, i, sparams.omeas,
                                                                true);
          }
        }
      }
    }

    void after_hmc_step(const size_t &i) {
      if (i > 0 && (i % sparams.N_save) == 0) { // saving U after each N_save trajectories
        std::string path_i = conf_path_basename + "." + std::to_string(i);
        if (sparams.do_mcmc) {
          U.save(path_i);
        }

        // online measurements
        bool do_omeas =
          sparams.do_omeas && i > sparams.omeas.icounter && i % sparams.omeas.nstep == 0;
        if (sparams.do_mcmc) {
          // check also if trajectory was accepted
          do_omeas = do_omeas && mdparams.getaccept();
        }

        if (do_omeas) {
          if (i == g_icounter&& sparams.do_mcmc) {
            return; // online measurements already done
          }

          /* loading the gauge configuration */
          int lerr = U.load(path_i);
          if (lerr == 1) {
            return;
          }
        }
        this->do_omeas_i(i, do_omeas);

        if (sparams.do_mcmc) { // storing last conf index (only after online measurements
                               // has been done)
          io::hmc::update_nconf_counter(sparams.conf_dir, g_heat, i, path_i);
        }
      }
    }

    void run(int ac, char *av[]) {
      this->print_program_info();
      this->print_git_info();

      this->parse_command_line(ac, av);
      this->parse_input_file();

      this->create_directories();

      this->init_gauge_conf();
      this->open_output_data();

      // Molecular Dynamics parameters
      md_params mdparams0(sparams.n_steps, sparams.tau);
      mdparams = mdparams0;

      this->init_monomials();

      // setting up the integrator
      md_integ = set_integrator<double, _u1>(sparams.integrator, sparams.exponent);

      // header: column names in the output
      std::string head_str = io::hmc::get_header(" ");
      std::cout << head_str;
      os << head_str;

      const std::string conf_path_basename = io::get_conf_path_basename(pparams, sparams);

      for (size_t i = g_icounter ; i < sparams.n_meas + g_icounter; i++) {
        this->do_hmc_step(i);
        this->after_hmc_step(i);
      }

      if (sparams.do_mcmc) {
        std::cout << "## Acceptance rate: " << rate / static_cast<double>(sparams.n_meas)
                  << std::endl;
        std::string path_final = conf_path_basename + ".final";
        U.save(path_final);
      }

      return;
    }
  };

} // namespace u1
