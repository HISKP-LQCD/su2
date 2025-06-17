/**
 * @file nested_sampling.hpp
 * @author Simone Romiti (simone.romiti.1994@gmail.com)
 * @brief class for Nested Sampling Algorithm. It works for integrals of functions of the
 * average plaquette only
 * @version 0.1
 * @date 2022-09-01
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdio.h>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xcsv.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xtensor.hpp>

#include "base_program.hpp"
#include "errors.hpp"
#include "io.hh"
#include "uniform_sweeps.hpp"

template <class Group>
class nested_sampling_algo : public base_program<Group, gp::nested_sampling> {
private:
  std::vector<size_t> indices; // list of configuration indices
  std::vector<double> Pi; // list of plaquette values for n_live points

  std::mt19937 engine; // engine for random number generation
  std::vector<std::mt19937> engines; // engines for exceptional overrelaxation updates

  std::string conf_counter_file; // file to save the configuration counter
  std::string output_data_file; // file to save the output.data
  std::ofstream os_nlive_conf; // configurations of the n_live points
  std::string path_nlive_conf; // path of os_nlive_conf
  std::ofstream os_nlive_idx; // configuration indices of the n_live points
  std::string path_nlive_idx; // path of os_nlive_idx

  int i_last_conf = 0; // index of the last configuration saved
  int i_step = 0; // index of the last NS step
  int i_dead = 0; // index of the dead point
  std::string conf_counter_path;
  std::string step_counter_path;

public:
  nested_sampling_algo() { (*this).algo_name = "nested_sampling"; }
  ~nested_sampling_algo() {
    os_nlive_conf.close();
    os_nlive_idx.close();
  }

  void print_program_info() const { std::cout << "## nested_sampling Algorithm\n"; }

  void save_nlive_status() {
    // saving the configuration of the n_live points
    io::vector_to_stream(Pi, (*this).os_nlive_conf, " ");
    io::vector_to_stream((*this).indices, (*this).os_nlive_idx, " ");
    return;
  }

  void parse_input_file(const YAML::Node &nd) {
    namespace in_nested_sampling = input_file_parsing::nested_sampling;
    in_nested_sampling::parse_input_file(nd, (*this).pparams, (*this).sparams);
    (*this).omeas = (*this).sparams.omeas;
    (*this).conf_path_basename =
      io::get_conf_path_basename((*this).pparams, (*this).sparams);
  }

  // unsorted list of plaquette values
  void init_nlive(const int &n_live, const int &seed) {
    (*this).Pi.resize(n_live);
    (*this).indices.resize(n_live);

    const double delta = 1.0; //(*this).sparams.delta;
    std::cout << "## Initial n_live values of the plaquette density (drawn at beta=0) \n";
    for (size_t i = 0; i < n_live; i++) {
      // creating a random gauge configuration
      hotstart<Group>((*this).U, seed + i, delta);
      std::string path_i = (*this).conf_path_basename + "." + std::to_string(i);
      (*this).U.save(path_i);

      // saving the value of the action (in this case, the plaquette density)
      const double pi = omeasurements::get_retr_plaquette_density((*this).U, "periodic");
      std::cout << pi << std::endl;
      Pi[i] = pi;
      (*this).indices[i] = i;
    }
    auto P_min_element = std::min_element(Pi.begin(), Pi.end());
    i_dead = std::distance(Pi.begin(), P_min_element); // conf. index of the minimum

    this->save_nlive_status();

    i_last_conf = n_live - 1; // index of the last configuration saved

    return;
  }

  std::vector<double> read_nlive_conf() {
    std::cout << "## Reading old n_live points from " << path_nlive_conf << std::endl;
    const int n_live = (*this).sparams.n_live;
    check_file_exists(path_nlive_conf, __func__);
    // std::ifstream in_file;
    // in_file.open(path_nlive_conf);
    check_file_exists(path_nlive_idx, __func__);
    // (*this).Pi.resize(n_live);
    // (*this).indices.resize(n_live);

    (*this).Pi = io::string_to_vector<double>(io::read_last_line(path_nlive_conf), " ");
    (*this).indices =
      io::string_to_vector<size_t>(io::read_last_line(path_nlive_idx), " ");

    // index of the last configuration saved
    i_last_conf = io::read_single_value<int>(conf_counter_file);

    return Pi;
  }

  void open_output_data() {
    conf_counter_file = (*this).sparams.conf_dir + "/conf_counter.txt";

    output_data_file =
      (*this).sparams.conf_dir + "/output." + (*this).algo_name + ".data";

    std::ios_base::openmode write_mode;
    if ((*this).sparams.continue_run == true) {
      write_mode = std::ios::app; // append to existing file
      check_file_exists(path_nlive_conf, __func__);
      check_file_exists(path_nlive_idx, __func__);
    } else {
      write_mode = std::ios::out; // create a new file
    }

    // index of configuration of dead point and value of P
    (*this).os.open(output_data_file, write_mode);
    (*this).os << std::scientific << std::setprecision(16);

    // list of n_live points, one for each NS step
    (*this).os_nlive_conf.open(path_nlive_conf, write_mode);
    (*this).os_nlive_conf << std::scientific << std::setprecision(16);

    // list of n_live configuration indices
    (*this).os_nlive_idx.open(path_nlive_idx, write_mode);
  }

  std::string get_path_conf(const int &i) const {
    return (*this).conf_path_basename + "." + std::to_string(i);
  }

  void do_omeas_i(const size_t &i) {
    namespace fsys = boost::filesystem;

    gaugeconfig<Group> U_i = (*this).U;

    if (!(*this).sparams.do_mcmc) { // doing only offline measurements
      const std::string path_i = get_path_conf(i);
      int ierrU = U_i.load(path_i);

      if (ierrU == 1) { // cannot load gauge config
        return; // simply ignore configuration
      }
    }

    if ((*this).omeas.polyakov.measure_it) {
      if ((*this).omeas.verbosity > 0) {
        std::cout << "## online measuring: Polyakov loop\n";
      }

      // creating the output directory
      std::ostringstream oss;
      oss << (*this).omeas.res_dir + "/" + (*this).omeas.polyakov.subdir << "/";
      std::string out_dir = oss.str();
      fsys::create_directories(fsys::absolute(out_dir)); // creating directory

      size_t i_orlx = 0;
      size_t N_overrelaxation = (*this).sparams.n_overrelaxation;
      const int n_threads = (*this).threads;
      while (i_orlx < N_overrelaxation + 1) {
        // initializing the engines
        for (size_t i_engine = 0; i_engine < n_threads; i_engine++) {
          (*this).engines[i_engine].seed(i * N_overrelaxation + i_orlx);
        }

        overrelaxation(U_i, (*this).engines, 1.0, false);

        std::string output_data_file = out_dir + "/Ploops." +
                                       boost::lexical_cast<std::string>(i) + "-" +
                                       boost::lexical_cast<std::string>(i_orlx);
        omeasurements::meas_polyakov(U_i, output_data_file);
        i_orlx++;
      }
    }

    return;
  }

  void offline_measurements() {
    const std::string output_data_file =
      (*this).sparams.conf_dir + "/output." + (*this).algo_name + ".data";

    std::ifstream file(output_data_file);

    if (!file.is_open()) {
      std::cerr << "Error opening: " << output_data_file << std::endl;
      std::abort();
    }

    std::string line;
    while (getline(file, line)) {
      std::istringstream iss(line);
      std::string firstColumn;
      if (getline(iss, firstColumn, ' ')) {
        const size_t i_conf = stod(firstColumn);

        this->do_omeas_i(i_conf);
      }
    }

    file.close();
  }

  void run(const YAML::Node &nd) {
    this->pre_run(nd); // prepare the algorithm

    path_nlive_conf = (*this).sparams.conf_dir + "/nlive_conf.dat";
    path_nlive_idx = (*this).sparams.conf_dir + "/nlive_idx.dat";
    // conf_counter_path = (*this).sparams.conf_dir + "/conf_counter.txt";
    // step_counter_path = (*this).sparams.conf_dir + "/step_counter.txt";

    bool do_omeas = (*this).sparams.do_omeas;
    bool do_mcmc = nd["nested_sampling"]["do_mcmc"].as<bool>();
    if (do_omeas && (!do_mcmc)) {
      this->offline_measurements();
      return; // do not run the algorithm, just measure observables
    }

    this->open_output_data(); // opening output files

    const size_t n_live = (*this).sparams.n_live;
    const size_t n_samples = (*this).sparams.n_samples;
    const size_t seed = (*this).sparams.seed;
    const double delta = (*this).sparams.delta;
    // number of sweeps per link, i.e. a multiple of the number of links
    const size_t n_sweeps_tot = ((*this).sparams.n_sweeps) * (*this).U.getSize();

    if ((*this).sparams.continue_run) {
      // i_step = 1 + this->read_from_counter(step_counter_path); // new step index
      read_nlive_conf();
    } else {
      std::cout << "## Initializing n_live points\n";
      init_nlive(n_live, seed);
    }

    if (do_omeas && !(*this).sparams.continue_run) {
      for (size_t j = 0; j < n_live; j++) {
        this->do_omeas_i((*this).indices[j]);
      }
    }

    // distribution of indices after the removal of one of the n_live points
    // ACHTUNG! right bound is included (it is the c++ syntax)
    std::uniform_int_distribution<> int_dist(0, n_live - 2);

    gaugeconfig<Group> &U_i = (*this).U; // configuration corresponding to that index

    // sampling n_samples points in the phase space
    for (size_t i = 0; i < n_samples; i++) {
      const int i_conf = i_last_conf + (1 + i); // configuration index

      // finding the minimum plaquette and appending it to the list
      auto P_min_element = std::min_element(Pi.begin(), Pi.end());
      i_dead = std::distance(Pi.begin(), P_min_element);
      const double Pmin = Pi[i_dead]; // minimum plaquette value

      // index of dead configuration
      const int i_dead_conf = (*this).indices[i_dead];
      ((*this).os) << i_dead_conf << " ";
      ((*this).os) << std::scientific << std::setprecision(16) << Pmin << std::endl;
      std::cout << i_dead_conf << " ";
      std::cout << std::scientific << std::setprecision(16) << Pmin << std::endl;

      // removing that element
      Pi.erase(Pi.begin() + i_dead);
      (*this).indices.erase((*this).indices.begin() + i_dead);

      // // Print indices and abort
      // std::cout << "i_dead_conf: " << i_dead_conf << std::endl;
      // std::cout << "i_last_conf: " << i_last_conf << std::endl;
      // std::cout << "i_conf: " << i_conf << std::endl;
      // std::cout << "Current indices: ";
      // for (const auto &idx : (*this).indices) {
      //   std::cout << idx << " ";
      // }
      // std::cout << std::endl;
      // std::abort();

      if ((*this).sparams.delete_dead_confs) {
        // removing dead configuration
        std::remove(this->get_path_conf(i_dead_conf).c_str());
      }

      // drawing a random element from the remained configurations
      std::mt19937 engine; // Random Number Generator (RNG)
      engine.seed(i_conf); // setting the seed of the RNG

      this->set_omp_threads();
      const int n_threads = (*this).threads;
      engines.resize(n_threads);

      const size_t ii_rand = int_dist(engine);
      const double Prand = Pi[ii_rand]; // value of the plaquette
      const size_t i_rand = (*this).indices[ii_rand]; // index of the configuration

      U_i.load(this->get_path_conf(i_rand), false, true);

      // applying a minimum of "n_sweeps_tot" sweeps to this configuration
      // to draw another one sampled from the constrained prior
      uniform_sweeps(U_i, Prand, Pmin, engine, delta, n_sweeps_tot);
      const double P_new = omeasurements::get_retr_plaquette_density(U_i, "periodic");

      // saving the new configuration of n_live points
      Pi.push_back(P_new);
      (*this).indices.push_back(i_conf);

      // saving the new configuration
      U_i.save(this->get_path_conf(i_conf));
      io::write_single_value<int>(i_conf, conf_counter_file);

      this->save_nlive_status(); // saving new configuration
      // this->write_to_counter(i_last_conf + i, conf_counter_path);
      // i_step++; // incrementing the step counter
      // this->write_to_counter(i_step, step_counter_path);

      if (do_omeas) {
        this->do_omeas_i(i_conf);
      }
    }

    (*this).os_nlive_conf.close();
    (*this).os_nlive_idx.close();
  }
};
