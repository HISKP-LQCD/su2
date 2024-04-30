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

#include <random>
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
      std::cout << pi << "\n";
      Pi[i] = pi;
      (*this).indices[i] = i;
    }
    std::cout << "---\n";
    return Pi;
  }

  std::vector<double> read_nlive_conf() {
    std::cout << "## Reading old n_live points from " << path_nlive_conf << "\n";
    const int n_live = (*this).sparams.n_live;
    check_file_exists(path_nlive_conf, __func__);
    std::vector<double> Pi(n_live);
    (*this).indices.resize(n_live);
    std::ifstream in_file;
    in_file.open(path_nlive_conf);
    const char delim = ' ';
    xt::xarray<double> conf_n_live = xt::load_csv<double>(in_file, delim, 1);

    std::cout << "Reading the following configuration\n";
    std::cout << conf_n_live << "\n";
    for (size_t i = 0; i < n_live; i++) {
      (*this).indices[i] = int(conf_n_live(i, 0));
      Pi[i] = conf_n_live(i, 1);

      // std::cout << i << " " << Pi[i] << "\n";
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
    } else {
      (*this).os.open(output_file, std::ios::out);
    }
    // opening the file for the last n_live points
    // check_file_exists(path_nlive_conf, __func__);
    (*this).os_nlive.open(path_nlive_conf, std::ios::out);
    // scientific notation's precision
    (*this).os << std::scientific << std::setprecision(16);
    (*this).os_nlive << std::scientific << std::setprecision(16);
  }

  void run(const YAML::Node &nd) {
    this->pre_run(nd);

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

    this->open_output_data();

    for (size_t i = 0; i < n_samples; i++) {
      // finding the minimum plaquette and appending it to the list
      const size_t i_min =
        std::distance(Pi.begin(), std::min_element(Pi.begin(), Pi.end()));

      const double Pmin = Pi[i_min];
      (*this).os << std::scientific << std::setprecision(16) << Pmin << "\n";
      // std::cout << i_last+i << " " << std::scientific << std::setprecision(16) << Pmin
      // << "\n";

      Pi.erase(Pi.begin() + i_min); // removing that element
      (*this).indices.erase((*this).indices.begin() + i_min);
      // drawing a random element from the remained configurations
      // random number generator
      std::mt19937 engine;
      engine.seed(i_last + i);
      std::uniform_int_distribution<> int_dist(0, (*this).indices.size());
      const size_t ii_rand = int_dist(engine);
      const double Prand = Pi[ii_rand]; // value of the plaquette
      const size_t i_rand = (*this).indices[ii_rand]; // index of the configuration
      gaugeconfig<Group> U_i = (*this).U; // configuration corresponding to that index
      U_i.load((*this).conf_path_basename + "." + std::to_string(i_rand), false, true);

      // applying n_sweeps_tot sweeps to this configuration
      // to draw another one sampled from the constrained prior
      uniform_sweeps(U_i, Prand, Pmin, engine, delta, n_sweeps_tot);
      const double P_new = omeasurements::get_retr_plaquette_density(U_i, "periodic");
      Pi.push_back(P_new);

      (*this).indices.push_back(i_last + i);
      // saving the new configuration
      std::string path_i = (*this).conf_path_basename + "." + std::to_string(i_last + i);
      U_i.save(path_i);
    }

    std::cout << "## Saving final configuration of n_live points\n";
    xt::xarray<double> conf_n_live(n_live);
    conf_n_live.resize({n_live, 2}); // i, P
    for (size_t i = 0; i < n_live; i++) {
      conf_n_live(i, 0) = (*this).indices[i];
      conf_n_live(i, 1) = Pi[i];
    }
    std::cout << conf_n_live << "\n";

    std::string header = "i P";
    io::xtensor_to_stream((*this).os_nlive, conf_n_live, " ", header);
    std::ofstream icounter((*this).sparams.conf_dir + "/icounter.txt");
    icounter << (i_last + n_samples);
    icounter.close();
  }
};
