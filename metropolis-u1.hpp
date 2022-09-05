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

  // namespace po = boost::program_options;
  // namespace gp = global_parameters;

  class metropolis_algo : public program<gp::metropolis_u1> {
private:

    std::vector<double> rate = {0., 0.};

  public:
    metropolis_algo() {}
    ~metropolis_algo() {}

    void print_program_info() const {
      std::cout << "## Metropolis Algorithm for U(1) gauge theory\n";
    }

    void parse_input_file() {
      namespace in_metropolis = input_file_parsing::u1::metropolis;
      int err = in_metropolis::parse_input_file(input_file, pparams, sparams);
      if (err > 0) {
        exit(1);
      }
      (*this).omeas = (*this).sparams.omeas;
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
      conf_path_basename = io::get_conf_path_basename(pparams, sparams);

      // load/set initial configuration
      const gaugeconfig<_u1> U0(pparams.Lx, pparams.Ly, pparams.Lz, pparams.Lt,
                                pparams.ndims, pparams.beta);
      U = U0;

      if (sparams.restart) {
        std::cout << "restart " << sparams.restart << std::endl;
        int err = U.load(conf_path_basename + "." + std::to_string(sparams.icounter));
        if (err != 0) {
          exit(1);
        }
      } else {
        hotstart(U, sparams.seed, sparams.heat);
      }

      // check gauge invariance, set up factors needed to normalise plaquette, spacial
      // plaquette
      double plaquette = this->gauge_energy<_u1>(pparams, U);

      const double fac = 1.0 / spacetime_lattice::num_pLloops_half(U.getndims());
      normalisation = fac / U.getVolume();
      facnorm = (pparams.ndims > 2) ? pparams.ndims / (pparams.ndims - 2) : 0;

      std::cout << "## Initial Plaquette: " << plaquette * normalisation << std::endl;

      random_gauge_trafo(U, 654321);

      // compute plaquette after the gauge transform
      plaquette = this->gauge_energy<_u1>(pparams, U);

      std::cout << "## Plaquette after rnd trafo: " << plaquette * normalisation
                << std::endl;
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

        double energy = this->gauge_energy<_u1>(pparams, U);

        double E = 0., Q = 0.;
        flat_spacetime::energy_density(U, E, Q);
        // measuring spatial plaquettes only means only (ndims-1)/ndims of all plaquettes
        // are measured, so need facnorm for normalization to 1
        std::cout << inew << " " << std::scientific << std::setw(18)
                  << std::setprecision(15) << energy * normalisation * facnorm << " ";
        os << inew << " " << std::scientific << std::setw(18) << std::setprecision(15)
           << energy * normalisation * facnorm << " ";

        energy = this->gauge_energy<_u1>(pparams, U);

        std::cout << energy * normalisation << " " << Q << " ";
        os << energy * normalisation << " " << Q << " ";
        this->energy_density(pparams, U, E, Q, false);
        std::cout << Q << std::endl;
        os << Q << std::endl;

        if (inew > 0 && (inew % sparams.N_save) == 0) {
          std::ostringstream oss_i;
          oss_i << conf_path_basename << "." << inew << std::ends;
          U.save(oss_i.str());
        }
      }

      return;
    }

    void save_acceptance_rates() {
      // save acceptance rates to additional file to keep track of measurements
      if (sparams.do_mcmc) {
        std::cout << "## Acceptance rate "
                  << rate[0] / double(sparams.n_meas)
                  << " temporal acceptance rate "
                  << rate[1] / double(sparams.n_meas) << std::endl;
        acceptancerates.open(sparams.conf_dir + "/acceptancerates.data", std::ios::app);
        acceptancerates << rate[0] / double(sparams.n_meas) << " "
                        << rate[1] / double(sparams.n_meas) << " "
                        << pparams.beta << " " << pparams.Lx << " " << pparams.Lt << " "
                        << pparams.xi << " " << sparams.delta << " " << sparams.heat
                        << " " << threads << " " << sparams.N_hit << " "
                        << sparams.n_meas << " " << sparams.seed << " " << std::endl;
        acceptancerates.close();

        std::ostringstream oss;
        oss << conf_path_basename << ".final" << std::ends;
        U.save(sparams.conf_dir + "/" + oss.str());
      }
      return;
    }

    /**
     * @brief do the i-th omeasurement
     *
     * @param i index of trajectory/sweep
     * @param inew
     */
    void do_omeas_i(const size_t &i, const size_t &inew) {
      if (sparams.do_omeas && inew != 0 && (inew % sparams.N_save) == 0) {
        if (!sparams.do_mcmc) {
          std::string path_i = conf_path_basename + "." + std::to_string(i);
          int ierrU = U.load(path_i);
          if (ierrU == 1) { // cannot load gauge config
            exit(1);
          }
        }
        if (omeas.potentialplanar || omeas.potentialnonplanar) {
          // smear lattice
          for (size_t smears = 0; smears < omeas.n_apesmear; smears += 1) {
            APEsmearing<double, _u1>(U, omeas.alpha, omeas.smear_spatial_only);
          }
          if (omeas.potentialplanar) {
            omeasurements::meas_loops_planar_pot(U, pparams, omeas.sizeWloops,
                                                 (*this).filename_coarse,
                                                 (*this).filename_fine, i);
          }

          if (omeas.potentialnonplanar) {
            omeasurements::meas_loops_nonplanar_pot(U, pparams, omeas.sizeWloops,
                                                    (*this).filename_nonplanar, i);
          }
        }
      }
      return;
    }



    void run(int ac, char *av[]) {
      this->print_program_info();
      this->print_git_info();

      this->parse_command_line(ac, av);
      this->parse_input_file();

      this->create_directories();

      this->set_omp_threads();
      this->set_potential_filenames();

      this->init_gauge_conf();
      this->open_output_data();

      /**
       * do measurements:
       * sweep: do N_hit Metropolis-Updates of every link in the lattice
       * calculate plaquette, spacial plaquette, energy density with and without cloverdef
       * and write to stdout and output-file save every nave configuration
       * */
      for (size_t i = sparams.icounter;
           i < sparams.n_meas * threads + sparams.icounter; i += threads) {
        // inew counts loops, loop-variable needed to have one RNG per thread with
        // different seeds for every measurement
        size_t inew = (i - sparams.icounter) / threads + sparams.icounter;

        this->do_sweep(i, inew);
        this->do_omeas_i(i, inew);
      }

      this->save_acceptance_rates();
    }
  };

} // namespace u1
