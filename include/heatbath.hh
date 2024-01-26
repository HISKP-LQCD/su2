/**
 * @file heatbath.hpp
 * @author Carsten Urbach (urbach@hiskp.uni-bonn.de)
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief
 * @version 0.1
 * @date 2022-05-30
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include "accum_type.hh"
#include "gaugeconfig.hh"
#include "get_staples.hh"
#include "random_element.hh"

#ifdef _USE_OMP_
#include <omp.h>
#endif

#include <algorithm>
#include <random>
#include <vector>

// function defined in the "proposed method" of
// - https://www.sciencedirect.com/science/article/pii/092056329290356W
// - https://arxiv.org/pdf/hep-lat/9210016v1.pdf
namespace hattori_nakajima {

  const double a_star = 0.798953686083986;
  const double epsilon = 0.001;

  // computes it from da = a-a^*, not from "a" itself
  double delta(const double &da) {
    return 0.35 * std::max(0.0, da) + 1.03 * sqrt(std::max(0.0, da));
  }

  double alpha(const double &a) {
    const double arg1 = sqrt(a * (2.0 - epsilon));
    const double arg2 = std::max(sqrt(epsilon * a), delta(a - a_star));
    return std::min(arg1, arg2);
  }

  double beta(const double &al, const double &a) {
    const double arg1 = std::pow(al, 2) / a;
    const double arg2 = (cosh(M_PI * al) - 1.0) / (exp(2 * a) - 1.0);
    return std::max(arg1, arg2) - 1.0;
  }

  double h(const double &al, const double &be, const double &x) {
    const double beta_fact = sqrt((1.0 - be) / (1.0 + be));
    const double A = (2.0 * x - 1.0) * atan(beta_fact * tanh(M_PI * al / 2.0));
    const double B = (1.0 / beta_fact) * tan(A);
    return (2.0 / al) * atanh(B);
  }

  double G(const double &al, const double &be, const double &a, const double &theta) {
    const double A1 = 1 - cos(theta);
    const double B1 = (cosh(al * theta) - 1.0) / (1.0 + be);
    return A1 - (1.0 / a) * log(1 + B1);
  }

  double g(const double &al, const double &be, const double &a, const double &x) {
    return exp(-a * G(al, be, a, h(al, be, x)));
  }

} // namespace hattori_nakajima

/**
 * @brief heatbath step
 *
 * Sampling each U_\mu(x) according to eq. (10) of
 * https://www.sciencedirect.com/science/article/pii/S0010465509002148   *
 *
 * @tparam URNG
 * @tparam Group
 * @param U gauge configuration
 * @param engine engine for RNG: size must be 2*number of threads
 * @param N_hit number of hits for the sampling
 * @param beta coupling beta in the action
 * @param xi bare anisotropy
 * @param anisotropic bool. flag, true when considering an anisotropic lattice
 * @return std::vector<double>: acceptance rates: {overall, only temporal ones}
 */
template <class URNG, class Group>
std::vector<double> heatbath(gaugeconfig<Group> &U,
                             std::vector<URNG> engine,
                             const double &beta,
                             const double &xi = 1.0,
                             const bool &anisotropic = false,
                             const bool &temporalonly = false, 
                             const bool &write_link = false);

template <class URNG>
std::vector<double> heatbath(gaugeconfig<u1> &U,
                             std::vector<URNG> engine,
                             const double &beta,
                             const double &xi = 1.0,
                             const bool &anisotropic = false,
                             const bool &temporalonly = false, 
                             const bool &write_link = false) {
  const double coupl_fact = (beta / double(U.getNc()));
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  typedef typename accum_type<u1>::type accum;
  double rate = 0.0, rate_time = 0.0, total_attempts = 0.0, temporal_attempts = 0.0;
  double changed_link_spatial = 0.0, changed_link_temporal = 0.0;

  const size_t endmu = temporalonly ? 1 : U.getndims(); 

  for (size_t x0_start = 0; x0_start < 2; x0_start++) {
#pragma omp parallel for reduction (+: rate, rate_time, total_attempts, temporal_attempts, changed_link_spatial, changed_link_temporal)
    for (size_t x0 = x0_start; x0 < U.getLt(); x0 += 2) {
      size_t thread_num = omp_get_thread_num();
      for (size_t x1 = 0; x1 < U.getLx(); x1++) {
        for (size_t x2 = 0; x2 < U.getLy(); x2++) {
          for (size_t x3 = 0; x3 < U.getLz(); x3++) {
            const std::vector<size_t> x = {x0, x1, x2, x3};
            for (size_t mu = 0; mu < endmu; mu++) {
              accum K;
              get_staples_MCMC_step(K, U, x, mu, xi, anisotropic);
              const double theta_stap = get_phase(K);
              const double rho = coupl_fact * get_abs(K);
              const double alpha = hattori_nakajima::alpha(rho);
              const double beta = hattori_nakajima::beta(alpha, rho);
              const double oldu = U(x, mu).geta();
              bool accept=false;
              double u1, u2;
              size_t attempt=0;
              while (!accept && attempt < 3000) {
                u1 = uniform(engine[thread_num]);
                u2 = uniform(engine[thread_num]);
                accept = (u2 < hattori_nakajima::g(alpha, beta, rho, u1));
                total_attempts+=1;
                if (mu == 0) {
                  temporal_attempts += 1;
                }
                attempt++;
              }
              if(attempt==3000) spacetime_lattice::fatal_error("too many attempts to generate new link!", __func__);
              U(x, mu).set(hattori_nakajima::h(alpha, beta, u1) - theta_stap);
              rate += 1;
              if (mu == 0) {
                rate_time += 1;
                changed_link_temporal += std::abs(U(x, mu).geta() - oldu);
              }
              else {
                changed_link_spatial += std::abs(U(x, mu).geta() - oldu);
              }
            }
          }
        }
      }
    }
  }
    if (write_link) {
      std::cout << std::setw(14) << "#### heatbath " << changed_link_spatial / double(U.getndims()*U.getVolume())
        << " " << changed_link_temporal / double(U.getVolume()) << std::endl;
    }

  const size_t den_acceptance_rate = temporalonly ? U.getVolume() : U.getSize();
  const std::vector<double> res = {double(rate) / double(den_acceptance_rate),
                                   double(rate_time) / double(U.getVolume()), 
                                   double(rate) / double(total_attempts), 
                                   double(rate_time) / double(temporal_attempts)};
  return res;
}

template <class URNG>
std::vector<double> heatbath_legacy(gaugeconfig<u1> &U,
                                    std::vector<URNG> engine,
                                    const double &beta,
                                    const double &xi = 1.0,
                                    const bool &anisotropic = false) {
  const double coupl_fact = (beta / double(U.getNc()));
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  typedef typename accum_type<u1>::type accum;
  double rate = 0.0, rate_time = 0.0;

  for (size_t x0_start = 0; x0_start < 2; x0_start++) {
#pragma omp parallel for
    for (size_t x0 = x0_start; x0 < U.getLt(); x0 += 2) {
      size_t thread_num = omp_get_thread_num();
      for (size_t x1 = 0; x1 < U.getLx(); x1++) {
        for (size_t x2 = 0; x2 < U.getLy(); x2++) {
          for (size_t x3 = 0; x3 < U.getLz(); x3++) {
            const std::vector<size_t> x = {x0, x1, x2, x3};
            for (size_t mu = 0; mu < U.getndims(); mu++) {
              accum K;
              get_staples_MCMC_step(K, U, x, mu, xi, anisotropic);
              const double theta_stap = get_phase(K);
              const double rho = coupl_fact * get_abs(K);
              const double phi1 = U(x, mu).geta() + theta_stap;
              const double phi2 = 2 * M_PI * uniform(engine[thread_num]);
              const double Sp = rho * (cos(phi2) - cos(phi1));
              bool accept = (Sp >= 0);
              if (!accept) {
                accept = (uniform(engine[thread_num]) < exp(Sp));
              }
              if (accept) {
                U(x, mu).set(phi2 - theta_stap);
                rate += 1;
                if (mu == 0) {
                  rate_time += 1;
                }
              }
            }
          }
        }
      }
    }
  }

  const std::vector<double> res = {double(rate) / double(U.getSize()),
                                   double(rate_time) / double(U.getVolume())};
  return res;
}

template <class URNG>
std::vector<double> heatbath(gaugeconfig<su2> &U,
                             std::vector<URNG> engine,
                             const double &beta,
                             const double &xi = 1.0,
                             const bool &anisotropic = false,
                             const bool &temporalonly = false, 
                             const bool &write_link = false) {
  spacetime_lattice::fatal_error("heatbath not implemented for SU(2)!", __func__);

  return {};
}

template <class URNG>
std::vector<double> heatbath(gaugeconfig<su3> &U,
                             std::vector<URNG> engine,
                             const double &beta,
                             const double &xi = 1.0,
                             const bool &anisotropic = false,
                             const bool &temporalonly = false, 
                             const bool &write_link = false) {
  spacetime_lattice::fatal_error("heatbath not implemented for SU(3)!", __func__);

  return {};
}
