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

#include "program-u1.hpp"

#include "detDDdag_monomial.hh"

namespace u1 {

  class measure_algo : public program<gp::measure_u1> {

  public:
    measure_algo() {}
    ~measure_algo() {}

    void print_program_info() const {
      std::cout << "## Measuring Tool for U(1) gauge theory" << std::endl;
    }

    void parse_input_file(const YAML::Node &nd) {
      namespace in_meas = input_file_parsing::u1::measure;
      in_meas::parse_input_file(nd, pparams, sparams);
      (*this).omeas = (*this).sparams;
      conf_path_basename = io::get_conf_path_basename(pparams, sparams);
    }


    void run(const YAML::Node &nd) {
      this->pre_run(nd);

      const size_t istart =
        omeas.icounter == 0 ? omeas.icounter + omeas.nstep : omeas.icounter;
      const size_t nmax = omeas.n_meas * omeas.nstep + omeas.icounter;
      for (size_t i = istart; i < nmax; i += omeas.nstep) {
        this->do_omeas_i(i);
      }

      return;
    }
  };

} // namespace u1
