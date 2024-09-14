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
#include "uniform_sweeps.hpp"

template <class Group>
class nested_sampling_algo : public base_program<Group, gp::nested_sampling> {
private:
  double rate = 0.0; // acceptance rate of the internal Metropolis update
  std::vector<size_t> indices; // list of gonfiguration indices
  // std::vector<double> plaquettes; // list of plaquette expectation values
  std::mt19937 engine; // engine for random number generation
  std::ofstream os_nlive; // output stram for the configuration of the final n_live points
  std::string path_nlive_conf; // path to configuration of the final n_live points
  int i_last; // index of the last configuration saved

public:
  nested_sampling_algo() { (*this).algo_name = "nested_sampling"; }
  ~nested_sampling_algo() { os_nlive.close(); }

  void print_program_info() const { std::cout << "## nested_sampling Algorithm\n"; }

  void parse_input_file(const YAML::Node &nd) {
    namespace in_nested_sampling = input_file_parsing::nested_sampling;
    in_nested_sampling::parse_input_file(nd, (*this).pparams, (*this).sparams);
    (*this).omeas = (*this).sparams.omeas;
    (*this).conf_path_basename =
      io::get_conf_path_basename((*this).pparams, (*this).sparams);
  }

  // unsorted list of plaquette values
  std::vector<double> init_nlive(const int &n_live, const int &seed) {
    double p0 = 0.0; // minimum value of plaquette --> minimum of e^(-S/beta)
    // finding the minimum value
    std::vector<double> Pi(n_live);
    (*this).indices.resize(n_live);

    const double delta = 1.0; //(*this).sparams.delta;
    std::cout << "## Initial n_live values of the plaquette\n";
    for (size_t i = 0; i < n_live; i++) {
      hotstart<Group>((*this).U, seed + i, delta);
      std::string path_i = (*this).conf_path_basename + "." + std::to_string(i);
      (*this).U.save(path_i);

      const double pi = omeasurements::get_retr_plaquette_density((*this).U, "periodic");
      std::cout << pi << std::endl;
      Pi[i] = pi;
      (*this).indices[i] = i;
    }
    std::cout << "---" << std::endl;
    return Pi;
  }

  std::vector<double> read_nlive_conf() {
    std::cout << "## Reading old n_live points from " << path_nlive_conf << std::endl;
    const int n_live = (*this).sparams.n_live;
    check_file_exists(path_nlive_conf, __func__);
    std::vector<double> Pi(n_live);
    (*this).indices.resize(n_live);
    std::ifstream in_file;
    in_file.open(path_nlive_conf);
    const char delim = ' ';
    xt::xarray<double> conf_n_live = xt::load_csv<double>(in_file, delim, 1);

    std::cout << "Reading the following configuration" << std::endl;
    std::cout << conf_n_live << std::endl;
    for (size_t i = 0; i < n_live; i++) {
      (*this).indices[i] = int(conf_n_live(i, 0));
      Pi[i] = conf_n_live(i, 1);
    }
    in_file.close();

    return Pi;
  }

  void open_output_data() {
    const std::string output_file =
      (*this).sparams.conf_dir + "/output." + (*this).algo_name + ".data";
    std::ofstream os_nlive;

    if ((*this).sparams.continue_run == true) {
      (*this).os.open(output_file, std::ios::app);
      // check if there is a saved configuration for the n_live points
      check_file_exists(path_nlive_conf, __func__);
    } else {
      (*this).os.open(output_file, std::ios::out);

      // opening the file for the last n_live points
      //(*this).os_nlive.open(path_nlive_conf, std::ios::out);
    }

    // scientific notation's precision
    (*this).os << std::scientific << std::setprecision(16);
    // (*this).os_nlive << std::scientific << std::setprecision(16);
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
      while (i_orlx < (*this).sparams.n_overrelaxation + 1) {
        overrelaxation(U_i, 1.0, false);

        std::string output_file = out_dir + "/Ploops." +
                                  boost::lexical_cast<std::string>(i) + "-" +
                                  boost::lexical_cast<std::string>(i_orlx);
        omeasurements::meas_polyakov(U_i, output_file);
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
    this->pre_run(nd);

    bool do_omeas = (*this).sparams.do_omeas;
    bool do_mcmc = nd["nested_sampling"]["do_mcmc"].as<bool>();
    if (!do_mcmc) {
      this->offline_measurements();
      return;
    }

    const size_t n_live = (*this).sparams.n_live;
    const size_t n_samples = (*this).sparams.n_samples;
    const size_t seed = (*this).sparams.seed;
    const double delta = (*this).sparams.delta;
    // number of sweeps over the entire lattice: e.g 1,2,3
    const size_t n_sweeps_tot = ((*this).sparams.n_sweeps) * (*this).U.getSize();

    path_nlive_conf = (*this).sparams.conf_dir + "/nlive_conf.data";
    std::vector<double> Pi;
    if ((*this).sparams.continue_run) {
      Pi = read_nlive_conf();

      // reading index of the last configuration
      std::ifstream icounter((*this).sparams.conf_dir + "/icounter.txt");
      icounter >> i_last; // Read the number from the file
      icounter.close(); // Close the input file

    } else {
      std::cout << "## Initializing n_live points\n";
      Pi = init_nlive(n_live, seed);
      i_last = n_live;
    }

    if (do_omeas && !(*this).sparams.continue_run) {
      for (size_t j = 0; j < n_live; j++) {
        this->do_omeas_i((*this).indices[j]);
      }
    }

    this->open_output_data();

    // distribution of indices after the removal of one of the n_live points
    // ACHTUNG! right bound is included (it is the c++ syntax)
    std::uniform_int_distribution<> int_dist(0, n_live - 2);

    gaugeconfig<Group> &U_i = (*this).U; // configuration corresponding to that index
    // sampling n_samples points in the phase space
    for (size_t i = 0; i < n_samples; i++) {
      const int i_conf = i_last + i; // configuration index
      // finding the minimum plaquette and appending it to the list
      auto P_min_element = std::min_element(Pi.begin(), Pi.end());
      const size_t i_min = std::distance(Pi.begin(), P_min_element);
      const double Pmin = Pi[i_min];

      // index of dead configuration
      const int i_dead_conf = (*this).indices[i_min];
      ((*this).os) << i_dead_conf << " ";
      ((*this).os) << std::scientific << std::setprecision(16) << Pmin << std::endl;
      std::cout << i_dead_conf << " ";
      std::cout << std::scientific << std::setprecision(16) << Pmin << std::endl;

      // removing that element
      Pi.erase(Pi.begin() + i_min);
      (*this).indices.erase((*this).indices.begin() + i_min);

      if ((*this).sparams.delete_dead_confs) {
        // removing dead configuration
        std::remove(this->get_path_conf(i_dead_conf).c_str());
      }

      // drawing a random element from the remained configurations
      std::mt19937 engine; // Random Number Generator (RNG)
      engine.seed(i_conf); // setting the seed of the RNG

      const size_t ii_rand = int_dist(engine);
      const double Prand = Pi[ii_rand]; // value of the plaquette
      const size_t i_rand = (*this).indices[ii_rand]; // index of the configuration

      U_i.load(this->get_path_conf(i_rand), false, true);

      // applying a minimum of "n_sweeps_tot" sweeps to this configuration
      // to draw another one sampled from the constrained prior
      uniform_sweeps(U_i, Prand, Pmin, engine, delta, n_sweeps_tot);
      const double P_new = omeasurements::get_retr_plaquette_density(U_i, "periodic");
      Pi.push_back(P_new);

      (*this).indices.push_back(i_conf);
      // saving the new configuration
      U_i.save(this->get_path_conf(i_conf));

      std::cout << "## Saving final configuration of n_live points\n";
      (*this).os_nlive.open(path_nlive_conf, std::ios::out);
      (*this).os_nlive << std::scientific << std::setprecision(16);
      (*this).os_nlive << "i P" << std::endl;
      for (size_t k = 0; k < n_live; k++) {
        (*this).os_nlive << (*this).indices[k] << " " << Pi[k] << std::endl;
      }
      (*this).os_nlive.close();

      std::ofstream icounter((*this).sparams.conf_dir + "/icounter.txt");
      icounter << (i_last + i);
      icounter.close();

      if (do_omeas) {
        this->do_omeas_i(i_conf);
      }
    }
  }
};