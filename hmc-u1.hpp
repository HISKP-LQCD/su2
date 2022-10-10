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
    double rate = 0.; // acceptance rate

    std::list<monomial<double, _u1> *> monomial_list; // list of monomials in the action

    flat_spacetime::gaugemonomial<double, _u1> *gm; // gauge monomial
    rotating_spacetime::gauge_monomial<double, _u1>
      *gm_rot; // gauge monomial with space rotation
    kineticmonomial<double, _u1> *km; // kinetic momomial (momenta)
    staggered::detDDdag_monomial<double, _u1> *detDDdag; // staggered fermions monomial

    // Molecular Dynamics (MD)
    integrator<double, _u1> *md_integ; // MD integrator
    md_params mdparams; // MD parameters

  public:
    hmc_algo() { algo_name = "hmc"; }
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

    void parse_input_file(const YAML::Node &nd) {
      namespace in_hmc = input_file_parsing::u1::hmc;
      in_hmc::parse_input_file(nd, pparams, sparams);
      (*this).omeas = (*this).sparams.omeas;
      conf_path_basename = io::get_conf_path_basename(pparams, sparams);
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

    void after_hmc_step(const size_t &i) {
      if (i > 0 && (i % sparams.N_save) == 0) { // saving U after each N_save trajectories
        std::string path_i = conf_path_basename + "." + std::to_string(i);
        if (sparams.do_mcmc) {
          U.save(path_i);
        }

        // online measurements
        bool do_omeas = sparams.do_omeas && (i > sparams.omeas.icounter) &&
                        ((i - sparams.omeas.icounter) <= sparams.omeas.n_meas) &&
                        (i % sparams.omeas.nstep == 0);
        if (sparams.do_mcmc) {
          // check also if trajectory was accepted
          do_omeas = do_omeas && mdparams.getaccept();
        }

        if (do_omeas) {
          this->do_omeas_i(i);
        }

        if (sparams.do_mcmc) { // storing last conf index (only after online measurements
                               // has been done)
          io::hmc::update_nconf_counter(sparams.conf_dir, g_heat, i, path_i);
        }
      }
    }

    void run(const YAML::Node &nd) {
      this->pre_run(nd);
      this->init_gauge_conf_mcmc();

      if (g_icounter == 0) {
        // header: column names in the output
        std::string head_str = io::hmc::get_header(" ");
        std::cout << head_str;
        os << head_str;
      }

      // Molecular Dynamics parameters
      md_params md_p0(sparams.n_steps, sparams.tau);
      mdparams = md_p0;

      this->init_monomials();

      // setting up the integrator
      md_integ = set_integrator<double, _u1>(sparams.integrator, sparams.exponent);

      for (size_t i = g_icounter; i < sparams.n_meas + g_icounter; i++) {
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
