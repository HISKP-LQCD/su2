/**
 * @file measure-u1.hpp
 * @author Simone Romiti (simone.romiti.1994@gmail.com)
 * @brief class for the offline measurements
 * @version 0.1
 * @date 2022-09-02
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "base_program.hpp"

#include "detDDdag_monomial.hh"

template <class Group> class measure_algo : public base_program<Group, gp::measure> {
public:
  measure_algo() {}
  ~measure_algo() {}

  void print_program_info() const {
    std::cout << "## Measuring Tool" << std::endl;
  }

  void parse_input_file(const YAML::Node &nd) {
    namespace in_meas = input_file_parsing::measure;
    in_meas::parse_input_file(nd, (*this).pparams, (*this).sparams);
    (*this).omeas = (*this).sparams;
    (*this).conf_path_basename =
      io::get_conf_path_basename((*this).pparams, (*this).sparams);
  }

  void run(const YAML::Node &nd) {
    this->pre_run(nd);

		std::string omeas_counter_file = (*this).omeas.res_dir + "/omeas_counter.txt";
    if ((*this).omeas.restart){
			std::ifstream omeas_counter(boost::filesystem::absolute(omeas_counter_file));
      omeas_counter >> (*this).omeas.icounter;
      omeas_counter.close();
    }

    size_t istart = (*this).omeas.icounter;
    this->set_potential_filenames();
    
    const size_t nmax = (*this).omeas.n_meas + (*this).omeas.icounter;
    for (size_t i = istart; i < nmax; i += (*this).omeas.nstep) {
      this->do_omeas_i(i);
			std::ofstream omeas_counter(omeas_counter_file);
      omeas_counter << i;
      omeas_counter.close();
    }

    return;
  }
};
