/**
 * @file metropolis.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief class for Metropolis Algorithm
 * @version 0.1
 * @date 2022-09-01
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "base_program.hpp"

template <class Group>
class metropolis_algo : public base_program<Group, gp::metropolis> {
private:
  std::vector<double> rate = {0., 0.};

public:
  metropolis_algo() { (*this).algo_name = "metropolis"; }
  ~metropolis_algo() {}

  void print_program_info() const {
    std::cout << "## Metropolis Algorithm for U(1) gauge theory\n";
  }

  void parse_input_file(const YAML::Node &nd) {
    namespace in_metropolis = input_file_parsing::metropolis;
    in_metropolis::parse_input_file(nd, (*this).pparams, (*this).sparams);
    (*this).omeas = (*this).sparams.omeas;
    (*this).conf_path_basename =
      io::get_conf_path_basename((*this).pparams, (*this).sparams);
  }

  void set_omp_threads() {
#ifdef _USE_OMP_
    /**
     * the parallelisation of the sweep-function first iterates over all odd points in t
     * and then over all even points because the nearest neighbours must not change
     * during the updates, this is not possible for an uneven number of points in T
     * */
    if ((*this).pparams.Lt % 2 != 0) {
      std::cerr << "For parallel computing an even number of points in T is needed!"
                << std::endl;
      omp_set_num_threads(1);
      std::cerr << "Continuing with one thread." << std::endl;
    }
    // set things up for parallel computing in sweep
    (*this).threads = omp_get_max_threads();
#else
    (*this).threads = 1;
#endif
    std::cout << "threads " << (*this).threads << std::endl;
  }

  template <class URNG>
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
      spacetime_lattice::fatal_error("Rotating metric not supported yet.", __func__);
      return {};
      // return rotating_spacetime::sweep(U, pparams.Omega, engines, delta, N_hit,
      //                                  pparams.beta, pparams.xi,
      //                                  pparams.anisotropic);
    } else {
      spacetime_lattice::fatal_error("Invalid metric.", __func__);
      return {};
    }
  }

  /**
   * @brief do the i-th sweep of the metropolis algorithm
   *
   * @param i trajectory index
   */
  void do_sweep(const size_t &i) {
    const size_t n_threads = (*this).threads;
    if ((*this).sparams.do_mcmc) {
      std::vector<std::mt19937> engines(n_threads);
      for (size_t i_engine = 0; i_engine < n_threads; i_engine++) {
        engines[i_engine].seed((*this).sparams.seed + i * n_threads + i_engine);
      }

      this->output_line(i);

      rate += this->sweep((*this).pparams, (*this).U, engines, (*this).sparams.delta,
                          (*this).sparams.N_hit, (*this).pparams.beta, (*this).pparams.xi,
                          (*this).pparams.anisotropic);

      if (i > 0 && (i % (*this).sparams.N_save) == 0) {
        std::ostringstream oss_i;
        oss_i << (*this).conf_path_basename << "." << i << std::ends;
        ((*this).U).save(oss_i.str());
      }
    }

    return;
  }

  // save acceptance rates to additional file to keep track of measurements
  void save_acceptance_rates() {
    const size_t n_threads = (*this).threads;
    if ((*this).sparams.do_mcmc) {
      std::cout << "## Acceptance rate " << rate[0] / double((*this).sparams.n_meas)
                << " temporal acceptance rate "
                << rate[1] / double((*this).sparams.n_meas) << std::endl;
      (*this).acceptancerates.open((*this).sparams.conf_dir + "/acceptancerates.data",
                                   std::ios::app);
      (*this).acceptancerates << rate[0] / double((*this).sparams.n_meas) << " "
                              << rate[1] / double((*this).sparams.n_meas) << " "
                              << (*this).pparams.beta << " " << (*this).pparams.Lx << " "
                              << (*this).pparams.Lt << " " << (*this).pparams.xi << " "
                              << (*this).sparams.delta << " " << (*this).sparams.heat
                              << " " << n_threads << " " << (*this).sparams.N_hit << " "
                              << (*this).sparams.n_meas << " " << (*this).sparams.seed
                              << " " << std::endl;
      (*this).acceptancerates.close();
    }
    return;
  }

  void run(const YAML::Node &nd) {
    this->pre_run(nd);
    this->init_gauge_conf_mcmc();
    this->open_output_data();
    this->set_omp_threads();
    if ((*this).sparams.do_omeas) {
      this->set_potential_filenames();
    }

    if ((*this).g_icounter == 0) {
      // header: column names in the output
      std::string head_str = io::get_header_1(" ");
      std::cout << head_str;
      (*this).os << head_str;
    }

    const size_t i_min = (*this).g_icounter;
    const size_t i_max = (*this).sparams.n_meas + (*this).g_icounter;
    /**
     * do measurements:
     * sweep: do N_hit Metropolis-Updates of every link in the lattice
     * calculate plaquette, spacial plaquette, energy density with and without cloverdef
     * and write to stdout and output-file save every nave configuration
     * */
    for (size_t i = i_min; i < i_max; i++) {
      this->do_sweep(i);
      bool do_omeas =
        ((*this).sparams.do_omeas && i != 0 && (i % (*this).sparams.N_save) == 0);
      this->after_MCMC_step(i, do_omeas);
    }
    this->save_acceptance_rates();
    this->save_final_conf();
  }
};
