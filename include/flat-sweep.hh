/**
 * @file flat-sweep.hpp
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

#include <random>
#include <vector>

namespace flat_spacetime {

  /**
   * @brief N_hit Metropolis-Updates
   * does N_hit Metropolis updates of every link (U -> R*U, where R is a random element):
   * - picks a point 'x' and direction '\mu'
   * - updates the link U_{\mu}(x) (not the surrounding ones)
   * - sums up the unchanged part in staples (be calculated once for every link) --> get
   * \Delta S
   * - picks the next point and direction, and repat until has gone through the entire
   * lattice repa
   *
   * Notes:
   * - For the update, the nearest neighbour links have to be constant, hence
   * parallelization is not trivial. It is done by first updating all even time slices and
   * then all odd time slices.
   * - The acceptance rate can be tuned with delta, which determines the possible regions
   * from which R is drawn.
   * - For the normalization of the temporal rate: There is only one temporal link for
   * each lattice point, so the normalization is done with U.getVolume() With this
   * definition, for xi=1 the acceptance rates only differ in the third significant digit,
   * so this is correct
   *
   * The change in the action \Delta S accepted with probability min(1, exp(-\Delta S)).
   * @tparam URNG
   * @tparam Group
   * @param U
   * @param engine
   * @param delta
   * @param N_hit
   * @param beta
   * @param xi bare anisotropy
   * @param anisotropic bool flag, true when considering an anisotropic lattice. In this
   * case the action weights the temporal (including links in direction 0) and spatial
   * links differently
   * @return std::vector<double> vector of links acceptance rate: {overall, only temporal ones}
   */
  template <class URNG, class Group>
  std::vector<double> sweep(gaugeconfig<Group> &U,
                            std::vector<URNG> engine,
                            const double &delta,
                            const size_t &N_hit,
                            const double &beta,
                            const double &xi = 1.0,
                            const bool &anisotropic = false) {
    std::uniform_real_distribution<double> uniform(0., 1.);
    typedef typename accum_type<Group>::type accum;
    size_t rate = 0, rate_time = 0;
#ifdef _USE_OMP_
#pragma omp parallel
    {
      size_t thread_num = omp_get_thread_num();
#else
    size_t thread_num = 0;
#endif
#pragma omp for reduction(+ : rate, rate_time)
      for (size_t x0 = 0; x0 < U.getLt(); x0 += 2) {
        // Cannot use elements of a vector as iteration variables in for-loop with OpenMP,
        // so use dummy variables
        /*
         * Is it more efficient to use 2 of every variable for measuring the rates, or
         * would it be better to use vectors for everything? When using vectors, a new
         * reduction directive would have to be declared for the vectors omp_priv has to
         * be initialized with the = operator, initialization to zero happens with
         * zerovector-function #pragma omp declare reduction (+ : std::vector<double> :
         * omp_out += omp_in) initializer (omp_priv = zerovector(omp_orig.size())) OpenMP
         * specifications say: If the initializer-expr is a function name with an argument
         * list, then one of the arguments must be omp_priv or the address of omp_priv
         * This is not the case here, but in an example the code compiled anyway.
         */

        Group R;
        for (size_t x1 = 0; x1 < U.getLx(); x1++) {
          for (size_t x2 = 0; x2 < U.getLy(); x2++) {
            for (size_t x3 = 0; x3 < U.getLz(); x3++) {
              std::vector<size_t> x = {x0, x1, x2, x3};
              for (size_t mu = 0; mu < U.getndims(); mu++) {
                accum K;
                get_staples(K, U, x, mu, xi, anisotropic);
                for (size_t n = 0; n < N_hit; n++) {
                  random_element(R, engine[thread_num], delta);
                  double deltaS = beta / static_cast<double>(U.getNc()) *
                                  (retrace(U(x, mu) * K) - retrace(U(x, mu) * R * K));
                  bool accept = (deltaS < 0);
                  if (!accept)
                    accept = (uniform(engine[thread_num]) < exp(-deltaS));
                  if (accept) {
                    U(x, mu) = U(x, mu) * R;
                    U(x, mu).restoreSU();
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
#pragma omp for reduction(+ : rate, rate_time)
      for (size_t x0 = 1; x0 < U.getLt(); x0 += 2) {
        Group R;
        for (size_t x1 = 0; x1 < U.getLx(); x1++) {
          for (size_t x2 = 0; x2 < U.getLy(); x2++) {
            for (size_t x3 = 0; x3 < U.getLz(); x3++) {
              std::vector<size_t> x = {x0, x1, x2, x3};
              for (size_t mu = 0; mu < U.getndims(); mu++) {
                accum K;
                get_staples(K, U, x, mu, xi, anisotropic);
                for (size_t n = 0; n < N_hit; n++) {
                  random_element(R, engine[thread_num], delta);
                  double deltaS = beta / static_cast<double>(U.getNc()) *
                                  (retrace(U(x, mu) * K) - retrace(U(x, mu) * R * K));
                  bool accept = (deltaS < 0);
                  if (!accept)
                    accept = (uniform(engine[thread_num]) < exp(-deltaS));
                  if (accept) {
                    U(x, mu) = U(x, mu) * R;
                    U(x, mu).restoreSU();
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
#ifdef _USE_OMP_
    }
#endif
    std::vector<double> res = {double(rate) / double(N_hit) / double(U.getSize()),
                               double(rate_time) / double(N_hit) / double(U.getVolume())};
    return res;
  }

  /**
   * same as sweep, but only one single rng is passed to the function
   * hypothesis: drawing a pseudorandom number is a bottleneck in parallelization if only
   * one rng-engine supplies numbers to all threads solve this bottleneck by introducing a
   * vector of rng-engines, so there is one engine for each thread this was done in the
   * standard sweep function, this function is only for testing purposes and should be
   * deleted when the testing is concluded
   * */
  template <class URNG, class Group>
  std::vector<double> sweepone(gaugeconfig<Group> &U,
                               URNG &engine,
                               const double delta,
                               const size_t N_hit,
                               const double beta,
                               const double xi = 1.0,
                               bool anisotropic = false) {
    std::uniform_real_distribution<double> uniform(0., 1.);
    typedef typename accum_type<Group>::type accum;
    size_t rate = 0, rate_time = 0;
#ifdef _USE_OMP_
#pragma omp parallel
    {
#endif
#pragma omp for reduction(+ : rate, rate_time)
      for (size_t x0 = 0; x0 < U.getLt(); x0 += 2) {
        // Cannot use elements of a vector as iteration variables in for-loop with OpenMP,
        // so use dummy variables
        Group R;
        for (size_t x1 = 0; x1 < U.getLx(); x1++) {
          for (size_t x2 = 0; x2 < U.getLy(); x2++) {
            for (size_t x3 = 0; x3 < U.getLz(); x3++) {
              std::vector<size_t> x = {x0, x1, x2, x3};
              for (size_t mu = 0; mu < U.getndims(); mu++) {
                accum K;
                get_staples(K, U, x, mu, xi, anisotropic);
                for (size_t n = 0; n < N_hit; n++) {
                  random_element(R, engine, delta);
                  double deltaS = beta / static_cast<double>(U.getNc()) *
                                  (retrace(U(x, mu) * K) - retrace(U(x, mu) * R * K));
                  bool accept = (deltaS < 0);
                  if (!accept)
                    accept = (uniform(engine) < exp(-deltaS));
                  if (accept) {
                    U(x, mu) = U(x, mu) * R;
                    U(x, mu).restoreSU();
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
#pragma omp for reduction(+ : rate, rate_time)
      for (size_t x0 = 1; x0 < U.getLt(); x0 += 2) {
        Group R;
        for (size_t x1 = 0; x1 < U.getLx(); x1++) {
          for (size_t x2 = 0; x2 < U.getLy(); x2++) {
            for (size_t x3 = 0; x3 < U.getLz(); x3++) {
              std::vector<size_t> x = {x0, x1, x2, x3};
              for (size_t mu = 0; mu < U.getndims(); mu++) {
                accum K;
                get_staples(K, U, x, mu, xi, anisotropic);
                for (size_t n = 0; n < N_hit; n++) {
                  random_element(R, engine, delta);
                  double deltaS = beta / static_cast<double>(U.getNc()) *
                                  (retrace(U(x, mu) * K) - retrace(U(x, mu) * R * K));
                  bool accept = (deltaS < 0);
                  if (!accept)
                    accept = (uniform(engine) < exp(-deltaS));
                  if (accept) {
                    U(x, mu) = U(x, mu) * R;
                    U(x, mu).restoreSU();
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
#ifdef _USE_OMP_
    }
#endif
    std::vector<double> res = {double(rate) / double(N_hit) / double(U.getSize()),
                               double(rate_time) / double(N_hit) / double(U.getVolume())};
    return res;
  }

} // namespace flat_spacetime
