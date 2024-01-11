/**
 * @file heatbath_overrelaxation.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief class for heatbath_overrelaxation Algorithm
 * @version 0.1
 * @date 2022-09-01
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "base_program.hpp"
#include "heatbath.hh"
#include "overrelaxation.hh"
#include "parse_input_file.hh"

template <class Group>
class heatbath_overrelaxation_algo
  : public base_program<Group, gp::heatbath_overrelaxation> {
public:
  heatbath_overrelaxation_algo() { (*this).algo_name = "heatbath_overrelaxation"; }
  ~heatbath_overrelaxation_algo() {}

  void print_program_info() const {
    std::cout << "## Heatbath+Overrelaxation Algorithm\n";
  }

  void parse_input_file(const YAML::Node &nd) {
    namespace in_heatbath_overrelaxation = input_file_parsing::heatbath_overrelaxation;
    in_heatbath_overrelaxation::parse_input_file(nd, (*this).pparams, (*this).sparams);
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
  void heatbath(const gp::physics &pparams,
                gaugeconfig<Group> &U,
                std::vector<URNG> engines,
                const size_t &N_hit,
                const double &beta,
                const double &xi = 1.0,
                const bool &anisotropic = false) {
    flat_spacetime::heatbath(U, engines, N_hit, pparams.beta, pparams.xi,
                             pparams.anisotropic);
  }

  void open_output_data() {
    const std::string out_data_file = "output.heatbath_overrelaxation.data";
    if ((*this).g_icounter == 0) {
      (*this).os.open((*this).sparams.conf_dir + "/" + out_data_file, std::ios::out);
    } else {
      (*this).os.open((*this).sparams.conf_dir + "/" + out_data_file, std::ios::app);
    }
  }

  /**
   * @brief do the i-th sweep of the heatbath_overrelaxation algorithm
   *
   * @param i trajectory index
   * @param inew
   */
  void do_heatbath(const size_t &i) {
    if ((*this).sparams.do_mcmc) {
      const int n_engines = (*this).threads;
      std::vector<std::mt19937> engines(n_engines);
      for (size_t i_engine = 0; i_engine < n_engines; i_engine += 1) {
        engines[i_engine].seed((*this).sparams.seed + i + i_engine);
      }

      this->heatbath((*this).pparams, (*this).U, engines, (*this).sparams.N_hit,
                     (*this).pparams.beta, (*this).pparams.xi,
                     (*this).pparams.anisotropic);
    }
  }

  void do_overrelaxation() {
    flat_spacetime::overrelaxation((*this).U, (*this).pparams.beta, (*this).pparams.xi,
                                   (*this).pparams.anisotropic);
  }

  void run(const YAML::Node &nd) {
    this->pre_run(nd);
    this->init_gauge_conf_mcmc();
    this->open_output_data();
    this->set_omp_threads();
    if ((*this).sparams.do_omeas) {
      this->set_potential_filenames();
    }

    (*this).os << "i E Q E_ss Q_ss\n";
    size_t i_min = (*this).g_icounter;
    size_t i_max = (*this).sparams.n_meas * ((*this).threads) + (*this).g_icounter;
    size_t i_step = (*this).threads; // avoids using the same RNG seed 
    /**
     * do measurements:
     * sweep: do N_hit heatbath_overrelaxation-Updates of every link in the lattice
     * calculate plaquette, spacial plaquette, energy density with and without cloverdef
     * and write to stdout and output-file save every nave configuration
     * */
    for (size_t i = i_min; i < i_max; i += i_step) {
      // inew counts loops, loop-variable needed to have one RNG per thread with
      // different seeds for every measurement
      size_t inew = (i - (*this).g_icounter) / (*this).threads + (*this).g_icounter;

      double E = 0., Q = 0.;
      std::cout << inew;
      (*this).os << inew;
      for (bool ss : {false, true}) {
        this->energy_density((*this).pparams, (*this).U, E, Q, false, ss);
        std::cout << " " << std::scientific << std::setprecision(15) << E << " " << Q;
        (*this).os << " " << std::scientific << std::setprecision(15) << E << " " << Q;
      }
      std::cout << "\n";
      (*this).os << "\n";

      this->do_heatbath(i);
      this->do_overrelaxation();

      if (inew > 0 && (inew % (*this).sparams.N_save) == 0) {
        std::ostringstream oss_i;
        oss_i << (*this).conf_path_basename << "." << inew << std::ends;
        ((*this).U).save(oss_i.str());
      }

      bool b1 = (*this).sparams.do_omeas;
      bool b2 = inew != 0;
      bool b3 = (inew % (*this).sparams.N_save) == 0;
      bool do_omeas = (b1 && b2 && b3);
      this->after_MCMC_step(inew, do_omeas);

      std::ostringstream oss;
      oss << (*this).conf_path_basename << ".final" << std::ends;
      ((*this).U).save((*this).sparams.conf_dir + "/" + oss.str());
    }
  }
};
