/**
 * @file metropolis-u1.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief class for Metropolis Algorithm for U(1) gauge theory
 * @version 0.1
 * @date 2022-09-01
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "base_program.hpp"

template <class Group> class metropolis_algo : public base_program<Group, gp::metropolis> {
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

  void open_output_data() {
    if ((*this).g_icounter == 0) {
      (*this).os.open((*this).sparams.conf_dir + "/output.u1-metropolis.data",
                      std::ios::out);
    } else {
      (*this).os.open((*this).sparams.conf_dir + "/output.u1-metropolis.data",
                      std::ios::app);
    }
  }

  /**
   * @brief do the i-th sweep of the metropolis algorithm
   *
   * @param i trajectory index
   * @param inew
   */
  void do_sweep(const size_t &i, const size_t &inew) {
    if ((*this).sparams.do_mcmc) {
      std::vector<std::mt19937> engines((*this).threads);
      for (size_t engine = 0; engine < (*this).threads; engine += 1) {
        engines[engine].seed((*this).sparams.seed + i + engine);
      }

      rate += this->sweep((*this).pparams, (*this).U, engines, (*this).sparams.delta,
                          (*this).sparams.N_hit, (*this).pparams.beta, (*this).pparams.xi,
                          (*this).pparams.anisotropic);

      double E = 0., Q = 0., energy, spatialnorm;
      // number of plaquettes is different for spatial-spatial and total
      double facnorm = ((*this).pparams.ndims > 2) ? (*this).pparams.ndims / ((*this).pparams.ndims - 2) : 0;
      // total number of plaquettes, factor 2 because we only sum up mu>nu
      double normalisation = 2.0 / (*this).pparams.ndims / ((*this).pparams.ndims - 1) / (*this).U.getVolume() / double((*this).U.getNc());
      std::cout << inew;
      (*this).os << inew;
      for (bool ss : {false, true}) {
        this->energy_density((*this).pparams, (*this).U, E, Q, false, ss);
        energy = this->gauge_energy((*this).pparams, (*this).U, ss);
        spatialnorm = ss ? facnorm : 1.0;
        std::cout << " " << std::scientific << std::setprecision(15) << energy*normalisation*spatialnorm << " " << E << " " << Q;
        (*this).os << " " << std::scientific << std::setprecision(15) << energy*normalisation*spatialnorm << " " << E << " " << Q;
      }
      std::cout << "\n";
      (*this).os << "\n";

      if (inew > 0 && (inew % (*this).sparams.N_save) == 0) {
        std::ostringstream oss_i;
        oss_i << (*this).conf_path_basename << "." << inew << std::ends;
        ((*this).U).save(oss_i.str());
      }
    }

    return;
  }

  // save acceptance rates to additional file to keep track of measurements
  void save_acceptance_rates() {
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
                              << " " << (*this).threads << " " << (*this).sparams.N_hit
                              << " " << (*this).sparams.n_meas << " "
                              << (*this).sparams.seed << " " << std::endl;
      (*this).acceptancerates.close();

      std::ostringstream oss;
      oss << (*this).conf_path_basename << ".final" << std::ends;
      ((*this).U).save((*this).sparams.conf_dir + "/" + oss.str());
    }
    return;
  }

  void run(const YAML::Node &nd) {
    this->pre_run(nd);
    this->init_gauge_conf_mcmc();
    this->set_omp_threads();

    (*this).os << "## i P E Q P_ss E_ss Q_ss\n";
    size_t i_min = (*this).g_icounter;
    size_t i_max = (*this).sparams.n_meas * ((*this).threads) + (*this).g_icounter;
    size_t i_step = (*this).threads;
    /**
     * do measurements:
     * sweep: do N_hit Metropolis-Updates of every link in the lattice
     * calculate plaquette, spacial plaquette, energy density with and without cloverdef
     * and write to stdout and output-file save every nave configuration
     * */
    for (size_t i = i_min; i < i_max; i += i_step) {
      // inew counts loops, loop-variable needed to have one RNG per thread with
      // different seeds for every measurement
      size_t inew = (i - (*this).g_icounter) / (*this).threads + (*this).g_icounter;

      this->do_sweep(i, inew);
      bool do_omeas =
        ((*this).sparams.do_omeas && inew != 0 && (inew % (*this).sparams.N_save) == 0);
      this->after_MCMC_step(inew, do_omeas);
    }

    this->save_acceptance_rates();
  }
};
