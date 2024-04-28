/**
 * @file nested_sampling.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief class for Nested Sampling Algorithm. It works for integrals of functions of the
 * average plaquette only
 * @version 0.1
 * @date 2022-09-01
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "base_program.hpp"
#include "errors.hpp"

template <class Group>
class nested_sampling_algo : public base_program<Group, gp::nested_sampling> {
private:
  double rate = 0.0; // acceptance rate of the internal Metropolis update
  std::vector<double> plaquettes; // list of plaquette expectation values
  std::vector<std::string> paths; // paths of saved configurations

public:
  nested_sampling_algo() {
    (*this).algo_name = "nested_sampling";
    for (size_t i = 0; i < (*this).sparams.n_live; i++) {
      std::string path_i = (*this).conf_path_basename + "." + std::to_string(i);
      paths.push_back(path_i);
    }
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

  void init_nlive(const int &n_live, const int &seed) {
    double p0 = 0.0; // maximum value of plaquette --> minumum of e^(-S/beta)
    // finding the maximum value out of a
    for (size_t i = 0; i < n_live; i++) {
      hotstart<Group>((*this).U, seed, true);
      std::string path_i = (*this).conf_path_basename + "." + std::to_string(i);
      (*this).U.save(path_i);
      const double pi =
        omeasurements::get_retr_plaquette_density((*this).U, (*this).pparams.bc);
      if (pi > p0) {
        p0 = pi;
      }
      (*this).plaquettes.push_back(pi);
    }
    // std::sort((*this).plaquettes.begin(), (*this).plaquettes.end(),
    //           std::greater<double>()); // sorting plaquettes
  }

  void run(const YAML::Node &nd) {
    this->pre_run(nd);
    this->init_gauge_conf_mcmc();
    const size_t n_live = (*this).sparams.n_live;
    const size_t seed = (*this).sparams.seed;
    this->init_nlive(n_live, seed);

    std::vector<double> Pi = (*this).plaquettes;
    for (size_t i = 0; i < (*this).sparams.n_live; i++) {
      std::sort(Pi.begin(), Pi.end(), std::greater<double>()); // sorting plaquettes
      (*this).plaquettes[i] = Pi[0]; // selecting the max plaquette
      Pi.erase(Pi.begin()); // removing that element
    }

    // this->init_gauge_conf_mcmc();
    this->open_output_data();
    this->init_nlive((*this).sparams.n_live, (*this).sparams.seed);
  }
};
