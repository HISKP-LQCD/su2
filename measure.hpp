/**
 * @file measure-u1.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief class for the offline measurements for the U(1) gauge theory
 * @version 0.1
 * @date 2022-09-02
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "program.hpp"

#include "detDDdag_monomial.hh"

template <class Group> class measure_algo : public program<Group, gp::measure_u1> {
public:
  measure_algo() {}
  ~measure_algo() {}

  void print_program_info() const {
    std::cout << "## Measuring Tool for U(1) gauge theory" << std::endl;
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

    const size_t istart = (*this).omeas.icounter == 0
                            ? (*this).omeas.icounter + (*this).omeas.nstep
                            : (*this).omeas.icounter;
    const size_t nmax =
      (*this).omeas.n_meas * (*this).omeas.nstep + (*this).omeas.icounter;
    for (size_t i = istart; i < nmax; i += (*this).omeas.nstep) {
      this->do_omeas_i(i);
    }

    return;
  }
};
