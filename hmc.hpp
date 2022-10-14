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
#include "program.hpp"

//#include "rotating-gaugemonomial.hpp"

#include "detDDdag_monomial.hh"

template <class Group> class hmc_algo : public program<Group, gp::hmc_u1> {
private:
  double rate = 0.; // acceptance rate

  std::list<monomial<double, Group> *> monomial_list; // list of monomials in the action

  flat_spacetime::gaugemonomial<double, Group> *gm; // gauge monomial
  // rotating_spacetime::gauge_monomial<double, Group>
  //   *gm_rot; // gauge monomial with space rotation
  kineticmonomial<double, Group> *km; // kinetic momomial (momenta)
  staggered::detDDdag_monomial<double, Group> *detDDdag; // staggered fermions monomial

  // Molecular Dynamics (MD)
  integrator<double, Group> *md_integ; // MD integrator
  md_params mdparams; // MD parameters

public:
  hmc_algo() { (*this).algo_name = "hmc"; }
  ~hmc_algo() {
    free(gm);
    // free(gm_rot);
    free(km);
    free(detDDdag);
    free(md_integ);
  }

  void print_program_info() const {
    std::cout << "## HMC Algorithm for U(1) gauge theory\n";
  }

  void parse_input_file(const YAML::Node &nd) {
    namespace in_hmc = input_file_parsing::u1::hmc;
    in_hmc::parse_input_file(nd, (*this).pparams, (*this).sparams);
    (*this).omeas = (*this).sparams.omeas;
    (*this).conf_path_basename =
      io::get_conf_path_basename((*this).pparams, (*this).sparams);
  }

  // generate list of monomials
  void init_monomials() {
    (*this).km = new kineticmonomial<double, Group>(0);
    (*this).km->setmdpassive();
    (*this).monomial_list.push_back(km);

    if ((*this).pparams.include_gauge) {
      if ((*this).pparams.rotating_frame) {
        spacetime_lattice::fatal_error("Rotating metric not supported yet.", __func__);
        // (*this).gm_rot =
        //   new rotating_spacetime::gauge_monomial<double, Group>(0,
        //   (*this).pparams.Omega);
        // (*this).monomial_list.push_back(gm_rot);
      } else {
        (*this).gm =
          new flat_spacetime::gaugemonomial<double, Group>(0, (*this).pparams.xi);
        (*this).monomial_list.push_back(gm);
      }
    }

    if ((*this).pparams.include_staggered_fermions) { // including S_F (fermionic) in
                                                      // the action
      (*this).detDDdag = new staggered::detDDdag_monomial<double, Group>(
        0, (*this).pparams.m0, (*this).sparams.solver, (*this).sparams.tolerance_cg,
        (*this).sparams.seed_pf, (*this).sparams.solver_verbosity);
      (*this).monomial_list.push_back(detDDdag);
    }
  }

  void do_hmc_step(const int &i) {
    if ((*this).sparams.do_mcmc) {
      (*this).mdparams.disablerevtest();

      if (i > 0 && (*this).sparams.N_rev != 0 && (i) % (*this).sparams.N_rev == 0) {
        (*this).mdparams.enablerevtest();
      }
      // PRNG engine
      std::mt19937 engine((*this).sparams.seed + i);
      // perform the MD update

      md_update((*this).U, engine, mdparams, monomial_list, *md_integ);

      const double energy = flat_spacetime::gauge_energy((*this).U);
      double E = 0., Q = 0.;
      flat_spacetime::energy_density((*this).U, E, Q);
      rate += mdparams.getaccept();

      std::cout << i << " " << (*this).mdparams.getaccept() << " " << std::scientific
                << std::setw(18) << std::setprecision(15)
                << energy * (*this).normalisation << " " << std::setw(15)
                << (*this).mdparams.getdeltaH() << " " << std::setw(15)
                << rate / static_cast<double>(i + 1) << " ";

      if ((*this).mdparams.getrevtest()) {
        std::cout << (*this).mdparams.getdeltadeltaH();
      } else {
        std::cout << "NA";
      }
      std::cout << " " << Q << std::endl;

      (*this).os << i << " " << (*this).mdparams.getaccept() << " " << std::scientific
                 << std::setw(18) << std::setprecision(15)
                 << energy * (*this).normalisation << " " << std::setw(15)
                 << (*this).mdparams.getdeltaH() << " " << std::setw(15)
                 << rate / static_cast<double>(i + 1) << " ";
      if ((*this).mdparams.getrevtest()) {
        (*this).os << (*this).mdparams.getdeltadeltaH();
      } else {
        (*this).os << "NA";
      }
      (*this).os << " " << Q << std::endl;
    }
  }

  void after_hmc_step(const size_t &i) {
    if (i > 0 && (i % (*this).sparams.N_save) ==
                   0) { // saving (*this).U after each N_save trajectories
      std::string path_i = (*this).conf_path_basename + "." + std::to_string(i);
      if ((*this).sparams.do_mcmc) {
        (*this).U.save(path_i);
      }

      // online measurements
      bool do_omeas =
        (*this).sparams.do_omeas && (i > (*this).sparams.omeas.icounter) &&
        ((i - (*this).sparams.omeas.icounter) <= (*this).sparams.omeas.n_meas) &&
        (i % (*this).sparams.omeas.nstep == 0);
      if ((*this).sparams.do_mcmc) {
        // check also if trajectory was accepted
        do_omeas = do_omeas && mdparams.getaccept();
      }

      if (do_omeas) {
        this->do_omeas_i(i);
      }

      if ((*this).sparams.do_mcmc) { // storing last conf index (only after online
                                     // measurements
        // has been done)
        io::hmc::update_nconf_counter((*this).sparams.conf_dir, (*this).g_heat, i,
                                      path_i);
      }
    }
  }

  void run(const YAML::Node &nd) {
    this->pre_run(nd);
    this->init_gauge_conf_mcmc();

    if ((*this).g_icounter == 0) {
      // header: column names in the output
      std::string head_str = io::hmc::get_header(" ");
      std::cout << head_str;
      (*this).os << head_str;
    }

    // Molecular Dynamics parameters
    md_params md_p0((*this).sparams.n_steps, (*this).sparams.tau);
    mdparams = md_p0;

    this->init_monomials();

    // setting up the integrator
    md_integ =
      set_integrator<double, Group>((*this).sparams.integrator, (*this).sparams.exponent);

    for (size_t i = (*this).g_icounter; i < (*this).sparams.n_meas + (*this).g_icounter;
         i++) {
      this->do_hmc_step(i);
      this->after_hmc_step(i);
    }

    if ((*this).sparams.do_mcmc) {
      std::cout << "## Acceptance rate: "
                << rate / static_cast<double>((*this).sparams.n_meas) << std::endl;
      std::string path_final = (*this).conf_path_basename + ".final";
      (*this).U.save(path_final);
    }

    return;
  }
};
