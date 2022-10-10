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

#include "program-u1.hpp"

namespace u1 {

  class metropolis_algo : public program<gp::metropolis_u1> {
  private:
    std::vector<double> rate = {0., 0.};

  public:
    metropolis_algo() { algo_name = "metropolis"; }
    ~metropolis_algo() {}

    void print_program_info() const {
      std::cout << "## Metropolis Algorithm for U(1) gauge theory\n";
    }

    void parse_input_file(const YAML::Node &nd) {
      namespace in_metropolis = input_file_parsing::u1::metropolis;
      in_metropolis::parse_input_file(nd, pparams, sparams);
      (*this).omeas = (*this).sparams.omeas;
      conf_path_basename = io::get_conf_path_basename(pparams, sparams);
    }

    void set_omp_threads() {
#ifdef _USE_OMP_
      /**
       * the parallelisation of the sweep-function first iterates over all odd points in t
       * and then over all even points because the nearest neighbours must not change
       * during the updates, this is not possible for an uneven number of points in T
       * */
      if (pparams.Lt % 2 != 0) {
        std::cerr << "For parallel computing an even number of points in T is needed!"
                  << std::endl;
        omp_set_num_threads(1);
        std::cerr << "Continuing with one thread." << std::endl;
      }
      // set things up for parallel computing in sweep
      threads = omp_get_max_threads();
#else
      threads = 1;
#endif
      std::cout << "threads " << threads << std::endl;
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

    void open_output_data() {
      if (sparams.icounter == 0) {
        os.open(sparams.conf_dir + "/output.u1-metropolis.data", std::ios::out);
      } else {
        os.open(sparams.conf_dir + "/output.u1-metropolis.data", std::ios::app);
      }
    }

    /**
     * @brief do the i-th sweep of the metropolis algorithm
     *
     * @param i trajectory index
     * @param inew
     */
    void do_sweep(const size_t &i, const size_t &inew) {
      if (sparams.do_mcmc) {
        std::vector<std::mt19937> engines(threads);
        for (size_t engine = 0; engine < threads; engine += 1) {
          engines[engine].seed(sparams.seed + i + engine);
        }

        rate += this->sweep(pparams, U, engines, sparams.delta, sparams.N_hit,
                            pparams.beta, pparams.xi, pparams.anisotropic);

        double E = 0., Q = 0.;
        std::cout << inew;
        os << inew;
        for (bool ss : {false, true}) {
          this->energy_density(pparams, U, E, Q, false, ss);
          std::cout << " " << std::scientific << std::setprecision(15) << E << " " << Q;
          os << " " << std::scientific << std::setprecision(15) << E << " " << Q;
        }
        std::cout << "\n";
        os << "\n";

        if (inew > 0 && (inew % sparams.N_save) == 0) {
          std::ostringstream oss_i;
          oss_i << conf_path_basename << "." << inew << std::ends;
          U.save(oss_i.str());
        }
      }

      return;
    }

    // save acceptance rates to additional file to keep track of measurements
    void save_acceptance_rates() {
      if (sparams.do_mcmc) {
        std::cout << "## Acceptance rate " << rate[0] / double(sparams.n_meas)
                  << " temporal acceptance rate " << rate[1] / double(sparams.n_meas)
                  << std::endl;
        acceptancerates.open(sparams.conf_dir + "/acceptancerates.data", std::ios::app);
        acceptancerates << rate[0] / double(sparams.n_meas) << " "
                        << rate[1] / double(sparams.n_meas) << " " << pparams.beta << " "
                        << pparams.Lx << " " << pparams.Lt << " " << pparams.xi << " "
                        << sparams.delta << " " << sparams.heat << " " << threads << " "
                        << sparams.N_hit << " " << sparams.n_meas << " " << sparams.seed
                        << " " << std::endl;
        acceptancerates.close();

        std::ostringstream oss;
        oss << conf_path_basename << ".final" << std::ends;
        U.save(sparams.conf_dir + "/" + oss.str());
      }
      return;
    }

    void run(const YAML::Node &nd) {
      this->pre_run(nd);
      this->init_gauge_conf_mcmc();
      this->set_omp_threads();

      os << "i E Q E_ss Q_ss\n";
      /**
       * do measurements:
       * sweep: do N_hit Metropolis-Updates of every link in the lattice
       * calculate plaquette, spacial plaquette, energy density with and without cloverdef
       * and write to stdout and output-file save every nave configuration
       * */
      for (size_t i = sparams.icounter; i < sparams.n_meas * threads + sparams.icounter;
           i += threads) {
        // inew counts loops, loop-variable needed to have one RNG per thread with
        // different seeds for every measurement
        size_t inew = (i - sparams.icounter) / threads + sparams.icounter;

        this->do_sweep(i, inew);
        if (sparams.do_omeas && inew != 0 && (inew % sparams.N_save) == 0) {
          this->do_omeas_i(i);
        }
      }

      this->save_acceptance_rates();
    }
  };

} // namespace u1
