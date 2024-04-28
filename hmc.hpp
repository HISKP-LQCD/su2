/**
 * @file hmc-u1.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief class for the hmc algorithm
 * @version 0.1
 * @date 2022-09-02
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "errors.hpp"
#include "base_program.hpp"
#include "md_update.hh"

#include "obc_gaugemonomial.hh"

// #include "rotating-gaugemonomial.hpp"

#include "detDDdag_monomial.hh"

template <class Group> class hmc_algo : public base_program<Group, gp::hmc> {
private:
  double rate = 0.0; // acceptance rate

  std::list<monomial<double, Group> *> monomial_list; // list of monomials in the action

  flat_spacetime::gaugemonomial<double, Group> *gm = nullptr; // gauge monomial

  // gauge monomial with other (i.e. not periodic)
  // boundary conditions, e.g. "spatial_open"
  obc::gaugemonomial<double, Group> *obc_gm = nullptr;

  // rotating_spacetime::gauge_monomial<double, Group>
  //   *gm_rot = nullptr; // gauge monomial with space rotation

  kineticmonomial<double, Group> *km = nullptr; // kinetic momomial (momenta)
  staggered::detDDdag_monomial<double, Group> *detDDdag = nullptr; // staggered fermions

  // Molecular Dynamics (MD)
  integrator<double, Group> *md_integ = nullptr; // MD integrator
  md_params mdparams; // MD parameters

public:
  hmc_algo() { (*this).algo_name = "hmc"; }
  ~hmc_algo() {
    delete gm;
    delete obc_gm;
    // delete gm_rot;
    delete km;
    delete detDDdag;
    delete md_integ;
  }

  void print_program_info() const { std::cout << "## HMC Algorithm\n"; }

  void parse_input_file(const YAML::Node &nd) {
    namespace in_hmc = input_file_parsing::hmc;
    in_hmc::parse_input_file(nd, (*this).pparams, (*this).sparams);
    (*this).omeas = (*this).sparams.omeas;
    (*this).conf_path_basename =
      io::get_conf_path_basename((*this).pparams, (*this).sparams);
  }

  // generate list of monomials
  /**
   * @brief initializing monomials and adding them to the list
   *
   */
  void init_monomials() {
    (*this).km = new kineticmonomial<double, Group>(0);
    (*this).km->setmdpassive();
    (*this).monomial_list.push_back(km);

    if ((*this).pparams.include_gauge) {
      if ((*this).pparams.rotating_frame) {
        fatal_error("Rotating metric not supported yet.", __func__);
        // (*this).gm_rot =
        //   new rotating_spacetime::gauge_monomial<double, Group>(0,
        //   (*this).pparams.Omega);
        // (*this).monomial_list.push_back(gm_rot);
      } else {
        if ((*this).pparams.bc != "periodic") {
          obc::weights w((*this).pparams.bc, (*this).pparams.Lx, (*this).pparams.Ly,
                         (*this).pparams.Lz, (*this).pparams.Lt, (*this).pparams.ndims);
          (*this).obc_gm =
            new obc::gaugemonomial<double, Group>(0, w, (*this).pparams.xi);
          (*this).monomial_list.push_back(obc_gm);
        } else {
          (*this).gm =
            new flat_spacetime::gaugemonomial<double, Group>(0, (*this).pparams.xi);
          (*this).monomial_list.push_back(gm);
        }
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

      // do the MD update
      md_update((*this).U, engine, mdparams, monomial_list, *md_integ);

      double E = 0., Q = 0.;
      flat_spacetime::energy_density((*this).U, E, Q);

      rate += mdparams.getaccept();

      std::cout << i << " " << (*this).mdparams.getaccept() << " " << std::scientific
                << std::setw(18) << std::setprecision(15) << E << " "
                << std::setw(15) << (*this).mdparams.getdeltaH() << " " << std::setw(15)
                << rate / static_cast<double>(i + 1) << " ";

      if ((*this).mdparams.getrevtest()) {
        std::cout << (*this).mdparams.getdeltadeltaH();
      } else {
        std::cout << "NA";
      }
      std::cout << " " << Q << std::endl;

      (*this).os << i << " " << (*this).mdparams.getaccept() << " " << std::scientific
                 << std::setw(18) << std::setprecision(15) << E << " "
                 << std::setw(15) << (*this).mdparams.getdeltaH() << " " << std::setw(15)
                 << rate / static_cast<double>(i + 1) << " ";
      if ((*this).mdparams.getrevtest()) {
        (*this).os << (*this).mdparams.getdeltadeltaH();
      } else {
        (*this).os << "NA";
      }
      (*this).os << " " << Q << std::endl;
    }
  }

  void run(const YAML::Node &nd) {
    this->pre_run(nd);
    this->init_gauge_conf_mcmc();
    this->open_output_data();

    if ((*this).sparams.do_omeas) {
      this->set_potential_filenames();
    }

    // Molecular Dynamics parameters
    md_params md_p0((*this).sparams.n_steps, (*this).sparams.tau);
    mdparams = md_p0;

    this->init_monomials();

    // setting up the integrator
    md_integ =
      set_integrator<double, Group>((*this).sparams.integrator, (*this).sparams.exponent);

    if ((*this).g_icounter == 0) {
      // header: column names in the output
      std::string head_str = io::hmc::get_header(" ");
      std::cout << head_str;
      (*this).os << head_str;
    }

    const size_t omeas_icounter = (*this).sparams.omeas.icounter;
    const size_t omeas_nmeas = (*this).sparams.omeas.n_meas;
    const size_t omeas_nstep = (*this).sparams.omeas.nstep;
    for (size_t i = (*this).g_icounter; i < (*this).sparams.n_meas + (*this).g_icounter;
         i++) {
      this->do_hmc_step(i);
      // online measurements
      const size_t i_off = i - omeas_icounter; // offset removed
      bool do_omeas = (*this).sparams.do_omeas && (i_off > 0);
      do_omeas = do_omeas && (i_off <= omeas_nmeas);
      do_omeas = do_omeas && ((i_off % omeas_nstep) == 0);
      if ((*this).sparams.do_mcmc) {
        // check also if trajectory was accepted
        do_omeas = do_omeas && mdparams.getaccept();
      }
      this->after_MCMC_step(i, do_omeas);
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
