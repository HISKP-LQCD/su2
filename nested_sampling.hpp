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

#include "base_program.hpp"
#include "errors.hpp"
#include "uniform_sweeps.hpp"

template <class Group>
class nested_sampling_algo : public base_program<Group, gp::nested_sampling> {
private:
  double rate = 0.0; // acceptance rate of the internal Metropolis update
  std::vector<size_t> indices; // list of gonfiguration indices
  std::vector<double> plaquettes; // list of plaquette expectation values
  std::mt19937 engine; // engine for random number generation

public:
  nested_sampling_algo() {
    (*this).algo_name = "nested_sampling";
    // for (size_t i = 0; i < (*this).sparams.n_live; i++) {
    //   std::string path_i = (*this).conf_path_basename + "." + std::to_string(i);
    //   paths.push_back(path_i);
    // }
  }
  ~nested_sampling_algo() {}

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
    double p0 = 0.0; // minimum value of plaquette --> minumum of e^(-S/beta)
    // finding the minimum value
    std::vector<double> Pi;
    const double delta = 1.0; //(*this).sparams.delta;
    std::cout << "## Initial n_live value of the plaquette\n";
    for (size_t i = 0; i < n_live; i++) {
      hotstart<Group>((*this).U, seed + i, delta);
      std::string path_i = (*this).conf_path_basename + "." + std::to_string(i);
      (*this).U.save(path_i);

      const double pi =
        omeasurements::get_retr_plaquette_density((*this).U, (*this).pparams.bc);
      std::cout << pi << "\n";
      Pi.push_back(pi);
      (*this).indices.push_back(i);
    }
    std::cout << "---\n";
    return Pi;
  }

  void run(const YAML::Node &nd) {
    this->pre_run(nd);
    this->init_gauge_conf_mcmc();
    this->open_output_data();

    const size_t n_live = (*this).sparams.n_live;
    const size_t n_samples = (*this).sparams.n_samples;
    const size_t seed = (*this).sparams.seed;
    const double delta = (*this).sparams.delta;
    const size_t n_sweeps_lattice = ((*this).sparams.n_sweeps) * (*this).U.getSize();

    std::vector<double> Pi = this->init_nlive(n_live, seed);

    for (size_t i = 0; i < n_samples; i++) {
      // finding the minimum plaquette and appending it to the list
      const size_t i_min =
        std::distance(Pi.begin(), std::min_element(Pi.begin(), Pi.end()));

      // std::vector<double> Pi_sorted = Pi;
      // std::sort(Pi_sorted.begin(), Pi_sorted.end());
      // std::cout << "n_live values\n";
      // for (size_t k = 0; k < Pi.size(); k++) {
      //   std::cout << i << " " << k << " " << Pi_sorted[k] << "\n";
      // }
      // std::cout << "---\n";

      const double Pmin = Pi[i_min];
      // std::cout << "Pmin: " << Pmin << "\n";
      (*this).plaquettes.push_back(Pmin);
      (*this).os << std::scientific << std::setprecision(16) << Pmin << "\n";
      std::cout << i << " " << std::scientific << std::setprecision(16) << Pmin << "\n";

      Pi.erase(Pi.begin() + i_min); // removing that element
      (*this).indices.erase((*this).indices.begin() + i_min);
      // drawing a random element from the remained configurations
      // random number generator
      std::mt19937 engine;
      engine.seed(i);
      std::uniform_int_distribution<> int_dist(0, (*this).indices.size());
      const size_t ii_rand = int_dist(engine);
      const double Prand = Pi[ii_rand]; // value of the plaquette
      // std::cout << "Prand : " << Prand << "\n";
      size_t i_rand = (*this).indices[ii_rand]; // index of the configuration
      gaugeconfig<Group> U_i = (*this).U; // configuration corresponding to that index
      U_i.load((*this).conf_path_basename + "." + std::to_string(i_rand), false);
      // const double Prand_check =
      //   omeasurements::get_retr_plaquette_density(U_i, "periodic");
      // std::cout << "Prand_check : " << Prand_check << "\n";
      // applying n_sweeps_lattice sweeps to this configuration
      // to draw another one sampled from the constrained prior
      uniform_sweeps(U_i, Prand, Pmin, engine, delta * double(n_live) / double(i),
                     n_sweeps_lattice);
      const double P_new = omeasurements::get_retr_plaquette_density(U_i, "periodic");
      Pi.push_back(P_new);

      // std::cout << i << " : " << Pi.size() << "\n";
      (*this).indices.push_back(n_live + i);
      // saving the new configuration
      std::string path_i = (*this).conf_path_basename + "." + std::to_string(n_live + i);
      U_i.save(path_i);
    }
  }
};
