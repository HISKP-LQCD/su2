/**
 * @file program-u1.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief mother class for all the classes that handle the programs runs
 * @version 0.1
 * @date 2022-09-01
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "flat-energy_density.hh"
#include "flat-gauge_energy.hpp"
#include "flat-sweep.hh" // flat spacetime
#include "gaugeconfig.hh"
#include "io.hh"
#include "omeasurements.hpp"
#include "parse_input_file.hh"
#include "random_gauge_trafo.hh"
#include "rotating-energy_density.hpp" // rotating spacetime
#include "rotating-gauge_energy.hpp" // rotating spacetime
#include "rotating-sweep.hpp" // rotating spacetime
#include "su2.hh"
#include "u1.hh"
#include "vectorfunctions.hh"
#include "version.hh"

#ifdef _USE_OMP_
#include <omp.h>
#endif

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

namespace u1 {

  namespace po = boost::program_options;
  namespace gp = global_parameters;

template<class sparam_type>
  class program {
  protected:
    gp::physics pparams; // physics parameters
    sparam_type sparams; // specific parameters to the given run
    gp::measure_u1 omeas; // omeasurements parameters
    std::string input_file; // yaml input file path
    size_t threads;

    // filename needed for saving results from potential and potentialsmall
    std::string filename_fine;
    std::string filename_coarse;
    std::string filename_nonplanar;

    std::string conf_path_basename; // basename for configurations
    gaugeconfig<_u1> U; // gauge configuration evolved through the algorithm
    double normalisation;
    size_t facnorm;

    std::ofstream os;
    std::ofstream acceptancerates;

  public:
    program() {}
    ~program() {}

    virtual void print_program_info() const = 0;

    void print_git_info() const {
      std::cout << "## GIT branch " << GIT_BRANCH << " on commit \n\n";
    }

    void print_info() const {
      this->print_program_info();
      this->print_git_info();
    }

    /**
     * @brief parsing the command line for the main program.
     * Parameters available: "help" and "file"
     * see implementation for details
     * @param ac argc from the standard main() function
     * @param av argv from the standard main() function
     */
    void parse_command_line(int ac, char *av[]) {
      po::options_description desc("Allowed options");
      desc.add_options()("help,h", "produce this help message")(
        "file,f", po::value<std::string>(&input_file)->default_value("NONE"),
        "yaml input file");

      po::variables_map vm;
      po::store(po::parse_command_line(ac, av, desc), vm);
      po::notify(vm);

      if (vm.count("help")) {
        std::cout << desc << "\n";
        exit(0);
      }
      return;
    }

    virtual void parse_input_file() = 0;

    /**
     * @brief create the necessary output directories
     */
    void create_directories() {
      namespace fsys = boost::filesystem;
      fsys::create_directories(fsys::absolute(sparams.conf_dir));
      fsys::create_directories(fsys::absolute(omeas.res_dir));
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
      filename_fine = io::measure::get_filename_fine(pparams, omeas);
      filename_coarse = io::measure::get_filename_coarse(pparams, omeas);
      filename_nonplanar = io::measure::get_filename_nonplanar(pparams, omeas);

      // write explanatory headers into result-files, also check if measuring routine is
      // implemented for given dimension
      if (omeas.potentialplanar) {
        io::measure::set_header_planar(pparams, omeas, filename_coarse, filename_fine);
      }
      if (omeas.potentialnonplanar) {
        io::measure::set_header_nonplanar(pparams, omeas, filename_nonplanar);
      }
      return;
    }

    template <class Group>
    double gauge_energy(const gp::physics &pparams, const gaugeconfig<Group> &U) {
      if (pparams.flat_metric) {
        return flat_spacetime::gauge_energy(U);
      }
      if (pparams.rotating_frame) {
        return rotating_spacetime::gauge_energy(U, pparams.Omega);
      } else {
        spacetime_lattice::fatal_error("Invalid metric when calling: ", __func__);
        return {};
      }
    }

    template <class T>
    void energy_density(const gp::physics &pparams,
                        const gaugeconfig<T> &U,
                        double &E,
                        double &Q,
                        bool cloverdef = true) {
      if (pparams.flat_metric) {
        flat_spacetime::energy_density(U, E, Q, false);
      }
      if (pparams.rotating_frame) {
        rotating_spacetime::energy_density(U, pparams.Omega, E, Q, false);
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

    virtual void run(int ac, char *av[]) = 0;
  };

} // namespace u1
