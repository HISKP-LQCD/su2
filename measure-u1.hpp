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
  private:
    std::string conf_path_basename;

    size_t g_icounter = 0; // 1st configuration(trajectory) to load from
    double normalisation;

  public:
    measure_algo() {}
    ~measure_algo() {}

    void print_program_info() const {
      std::cout << "## Measuring Tool for U(1) gauge theory" << std::endl;
    }

    void parse_input_file() {
      namespace in_meas = input_file_parsing::u1::measure;
      in_meas::parse_input_file(input_file, pparams, sparams);
      (*this).omeas = (*this).sparams;
    }

    void init_gauge_conf() {
      gaugeconfig<_u1> U0(pparams.Lx, pparams.Ly, pparams.Lz, pparams.Lt, pparams.ndims,
                          pparams.beta);
      U = U0;
    }

    void run(int ac, char *av[]) {
      this->pre_run(ac, av);

      const size_t istart =
        omeas.icounter == 0 ? omeas.icounter + omeas.nstep : omeas.icounter;
      for (size_t i = istart; i < omeas.n_meas * omeas.nstep + omeas.icounter;
           i += omeas.nstep) {
        this->do_omeas_i(i);
      }

      return;
    }
  };

} // namespace u1
