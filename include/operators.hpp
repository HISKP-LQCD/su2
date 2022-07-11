// operators.hh
/**
 * @brief different operators on the lattice which can be measured for the use in
 * calculation of the glueball mass In contrast to the gauge_energy or Wilson-loop
 * operators, these operators should only ever look at one timeslice
 */

#pragma once

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "flat-gradient_flow.hh"
#include "gaugeconfig.hh"
#include "loops.hpp"
#include "parameters.hh"
#include "propagator.hpp"
#include "vectorfunctions.hh"
#include "wilsonloop.hh"

using Complex = std::complex<double>;
#pragma omp declare reduction (+ : std::complex<double> : omp_out += omp_in) initializer (omp_priv = 0)

namespace operators {

  template <class T> using nd_max_arr = spacetime_lattice::nd_max_arr<T>;

  /**
   * idea: make operators into a struct, so it is more compact. Maybe extend it into a
   * class with methods like checking for a closed loop and rotating an operator?
   * **/
  //~ struct operatorpath{
  //~ std::vector<size_t> lengths;
  //~ std::vector<size_t> directions;
  //~ std::vector<bool> sign;
  //~ } //struct operatorpath

  /**
   * @brief computes the sum of the trace of all plaquettes in the mu-nu-plane in the time
   * slice t
   * **/
  template <class Group>
  Complex plaquette_one_timeslice(const gaugeconfig<Group> &U,
                                  const size_t &t,
                                  const size_t &mu,
                                  const size_t &nu) {
    Complex res = 0.;

#pragma omp parallel for reduction(+ : res)
    for (size_t x = 0; x < U.getLx(); x++) {
      for (size_t y = 0; y < U.getLy(); y++) {
        for (size_t z = 0; z < U.getLz(); z++) {
          std::vector<size_t> vecx = {t, x, y, z};
          std::vector<size_t> vecxplusmu = vecx;
          std::vector<size_t> vecxplusnu = vecx;
          vecxplusmu[mu] += 1;
          vecxplusnu[nu] += 1;
          res += trace(U(vecx, mu) * U(vecxplusmu, nu) * U(vecxplusnu, mu).dagger() *
                       U(vecx, nu).dagger());
          vecxplusmu[mu] -= 1;
          vecxplusnu[nu] -= 1;
        }
      }
    }
    return res;
  }

  /**
   * @brief This multiplies all the link matrices along a given path and stores the
   * results in the group accum type The path starts at x, and is described by the vectors
   * lengths, directions, and sign
   * @param U configuration upon which the operator is measured
   * @param x vector of the starting point of the operator
   * @param lengths, directions, sign vectors describing the path of the operator
   * @param P parity, 0 for no parity operation, +/- for positive and negative parity
   * @note The vectors must have the same length. For each element, the path goes
   * lengths[i] steps into direction[i], with sign[i] determining if the path is traced
   * forwards or backwards. It makes sense to have all elements of direction in (0,
   * U.getndims-1), but if this is not the case, the internal function getindex will
   * ensure that directions[i] is interpreted as directions[i]%U.getndims
   * @return accum type of the group of the ordered product of all matrices in the
   * operator, so Complex for U(1) and SU(2) for SU(2). This makes it easier to take
   * either the real or imaginary part, with or without the trace Maybe it makes more
   * sense to return the trace in a complex number? Are there any cases where the trace is
   * not used in an operator?
   * @note example: P_mu,nu(x): lengths={1,1,1,1}, directions={mu, nu, mu, nu},
   * sign={true, true, false, false} example: W(t=1, x=4): lengths={1,4,1,4},
   * directions={0,1,0,1}, sign={true, true, false, false} example: chair xz-yz:
   * lengths={1,1,1,1,1,1}, directions={1,2,3,2,1,3}, sign={t,t,t,f,f,f} example: chair
   * xy-yz: lengths={1,1,1,1,1,1}, directions={1,3,2,3,2,1}, sign={t,t,t,f,f,f}
   * **/
  template <class Group>
  typename accum_type<Group>::type
  arbitrary_operator(const gaugeconfig<Group> &U,
                     const std::vector<size_t> &x,
                     const std::vector<size_t> &lengths,
                     const std::vector<size_t> &directions,
                     const std::vector<bool> &sign,
                     const size_t P = 0) {
    typedef typename accum_type<Group>::type accum;
    //~ if(P!=0){
    //~ std::cerr << "Parity conjugation is not yet implemented correctly! Your results
    // would be wrong. Aborting" << std::endl; ~ abort();
    //~ }
    if ((lengths.size() != directions.size()) || (lengths.size() != sign.size())) {
      std::cerr << "the lengths of the descriptors of the loop are not equal, no loop "
                   "can be calculated!"
                << std::endl;
      abort();
    }
    std::vector<size_t> xrun = x;
    accum L(1., 0.);
    for (size_t i = 0; i < lengths.size(); i++) {
      if (sign[i]) {
        for (size_t j = 0; j < lengths[i]; j++) {
          if (P != 0 && directions[i] > 0) {
            xrun[directions[i]]--;
            L *= U(xminusmu(invertspace(xrun), directions[i]), directions[i]).dagger();
          } else {
            if (P != 0) {
              L *= U(invertspace(xrun), directions[i]);
            } else {
              L *= U(xrun, directions[i]);
            }
            xrun[directions[i]]++;
          }
        }
      }
      if (!sign[i]) {
        for (size_t j = 0; j < lengths[i]; j++) {
          if (P != 0 && directions[i] > 0) {
            L *= U(xminusmu(invertspace(xrun), directions[i]), directions[i]);
            xrun[directions[i]]++;
          } else {
            xrun[directions[i]]--;
            if (P != 0) {
              L *= U(invertspace(xrun), directions[i]).dagger();
            } else {
              L *= U(xrun, directions[i]).dagger();
            }
          }
        }
      }
    }
    if (xrun != x) {
      std::cerr << "The loop was not closed!" << std::endl;
      //~ abort();
    }
    if (P < 0) {
      L = -L;
    }
    return L;
  }

  /**
   * @brief returns the sum of the real trace of the arbitrary operator
   * given by the vectors lengths, directions and sign in the timeslice with index t
   * see documentation of arbitrary_operator
   * **/
  template <class Group>
  double measure_re_arbitrary_loop_one_timeslice(const gaugeconfig<Group> &U,
                                                 const size_t t,
                                                 const std::vector<size_t> &lengths,
                                                 const std::vector<size_t> &directions,
                                                 const std::vector<bool> &sign) {
    double res = 0.;
    typedef typename accum_type<Group>::type accum;

#pragma omp parallel for reduction(+ : res)
    for (size_t x = 0; x < U.getLx(); x++) {
      for (size_t y = 0; y < U.getLy(); y++) {
        for (size_t z = 0; z < U.getLz(); z++) {
          std::vector<size_t> vecx = {t, x, y, z};
          res += retrace(arbitrary_operator(U, vecx, lengths, directions, sign));
        }
      }
    }
    return res;
  }

  /**
   * @brief returns the sum of the real trace of the arbitrary operator
   * given by the vectors lengths, directions and sign calculated in the entire lattice
   * see documentation of arbitrary_operator
   * **/
  template <class Group>
  double measure_re_arbitrary_loop_lattice(const gaugeconfig<Group> &U,
                                           const std::vector<size_t> &lengths,
                                           const std::vector<size_t> &directions,
                                           const std::vector<bool> &sign) {
    double res = 0.;
    typedef typename accum_type<Group>::type accum;

#pragma omp parallel for reduction(+ : res)
    for (size_t t = 0; t < U.getLt(); t++) {
      for (size_t x = 0; x < U.getLx(); x++) {
        for (size_t y = 0; y < U.getLy(); y++) {
          for (size_t z = 0; z < U.getLz(); z++) {
            std::vector<size_t> vecx = {t, x, y, z};
            res += retrace(arbitrary_operator(U, vecx, lengths, directions, sign));
          }
        }
      }
    }
    return res;
  }

  /**
   * @brief according to PC-specifications, returns an operator of the given loop with
   * parity, charge conjugation eigenvalues as given
   * **/
  template <class Group>
  double measure_arbitrary_loop_one_timeslice_PC(const gaugeconfig<Group> &U,
                                                 const size_t t,
                                                 const std::vector<size_t> &lengths,
                                                 const std::vector<size_t> &directions,
                                                 const std::vector<bool> &sign,
                                                 const bool P,
                                                 const bool C) {
    double res = 0.;
    typedef typename accum_type<Group>::type accum;
    accum K1, K2;
    size_t par = (P ? +1 : -1);
    //~ gaugeconfig<Group> PU=parityinvert(U);

#pragma omp parallel for reduction(+ : res)
    for (size_t x = 0; x < U.getLx(); x++) {
      for (size_t y = 0; y < U.getLy(); y++) {
        for (size_t z = 0; z < U.getLz(); z++) {
          std::vector<size_t> vecx = {t, x, y, z};
          //~ res += retrace(arbitrary_operator(U, vecx, lengths, directions, sign));
          K1 = arbitrary_operator(U, vecx, lengths, directions, sign, /*P=*/0);
          K2 = arbitrary_operator(U, vecx, lengths, directions, sign, /*P=*/par);
          if (C) { // C=+
            res += retrace(K1 + K2);
          } else { // C=-
            res += imtrace(K1 + K2);
          }
        }
      }
    }
    return res;
  }

  template <class T, class Group>
  T plaquette_Pij(const gaugeconfig<Group> &U,
                  const nd_max_arr<int> &x,
                  const size_t &i,
                  const size_t &j,
                  const bool &P) {
    T res;

    if (!P) {
      nd_max_arr<int> xpi = x;
      nd_max_arr<int> xpj = x;

      xpi[i]++;
      xpj[j]++;

      res = U(x, i) * U(xpi, j) * U(xpj, i).dagger() * U(x, j).dagger();

    } else {
      const nd_max_arr<int> Px = {x[0], -x[1], -x[2],
                                  -x[3]}; // spatial parity applied to 'x'
      nd_max_arr<int> Px_mimj = Px;
      nd_max_arr<int> Px_mi = Px;
      nd_max_arr<int> Px_mj = Px;

      Px_mimj[i]--;
      Px_mimj[j]--;
      Px_mi[i]--;
      Px_mj[j]--;

      res = U(Px_mi, i).dagger() * U(Px_mimj, j).dagger() * U(Px_mimj, i) * U(Px_mj, j);
    }
    return res / double(U.getNc());
  }

  /**
   * @brief \vec{p}=\vec{0} (at rest) trace of the spatial plaquette P*(U_ij) [P=spatial
   * parity operator]
   */
  template <class T, class Group>
  T rest_plaquette_P_ij(const gaugeconfig<Group> &U,
                        const size_t &t,
                        const size_t &i,
                        const size_t &j,
                        const bool &P) {
    T res;

#pragma omp parallel for reduction(+ : res)
    for (int x1 = 0; x1 < U.getLx(); x1++) {
      for (int x2 = 0; x2 < U.getLy(); x2++) {
        for (int x3 = 0; x3 < U.getLz(); x3++) {
          const nd_max_arr<int> x = {int(t), x1, x2, x3};
          res += plaquette_Pij<T, Group>(U, x, i, j, P);
        }
      }
    }
    return res / double(U.getNc()) / (double(U.getVolume()) / double(U.getLt()));
  }



  /**
   * @brief trace of product of links along a path, projected to '0' spatial momentum
   * 
   * @tparam T type to b returbed, e.g. std::complex<double>
   * @tparam Group 
   * @param t time
   * @param U gauge configuration
   * @param path_lat path along the lattice
   * @param P apply parity transformation
   * @return T trace of the product of the links along the path specified by variable 'path'
   */
  template <class T, class Group>
  T trace_rest_loop(const int& t, const gaugeconfig<Group> &U, const links::path &path_lat, const bool& P) {
    T res;
#pragma omp parallel for reduction(+ : res)
    for (int x1 = 0; x1 < U.getLx(); x1++) {
      for (int x2 = 0; x2 < U.getLy(); x2++) {
        for (int x3 = 0; x3 < U.getLz(); x3++) {
          const nd_max_arr<int> x = {t, x1, x2, x3};
          res += trace(path_lat.to_gauge<Group>(U, x, P));
        }
      }
    }
    return res / double(U.getNc()) / (double(U.getVolume()) / double(U.getLt()));
  }

  /**
   * @brief Get the sum tr U ij object
   * returns \frac{2}{(d-1)*(d-2)*Nc} \sum_{i<j} Tr(P*U_{ij}(t,
   * \vec{x}))
   */
  template <class T, class Group>
  T get_tr_sum_U_ij(const gaugeconfig<Group> &U,
                    const nd_max_arr<int> &x,
                    const bool &P) {
    T Uij;
#pragma omp parallel for reduction(+ : Uij) collapse(2)
    for (size_t i = 1; i < U.getndims(); i++) {
      for (size_t j = 2; j < U.getndims(); j++) {
        if (j <= i) { // we want j>i. If j<=i --> do nothing
          continue;
        }
        Uij += operators::plaquette_Pij<T, Group>(U, x, i, j, P);
      }
    }
    const double dims_fact = spacetime_lattice::num_pLloops_half(U.getndims() - 1);
    Uij /= double(dims_fact);
    return Uij;
  }

  /**
   * @brief Get the U ij object at rest
   * returns \frac{2}{(d-1)*(d-2) * Lx*Ly*Lz * Nc} \sum_{\vec{x}} \sum_{i<j}
   * Tr(P*U_{ij}(t, \vec{x})) P=parity operator
   * @tparam Group
   * @param U gauge configuration
   * @param P whether to apply the parity operator or not
   */
  template <class T, class Group>
  T get_rest_tr_sum_U_ij(const gaugeconfig<Group> &U, const size_t &t, const bool &P) {
    T Uij;
    for (size_t i = 1; i < U.getndims(); i++) {
      for (size_t j = i + 1; j < U.getndims(); j++) {
        Uij += operators::rest_plaquette_P_ij<T, Group>(U, t, i, j, P);
      }
    }
    const double dims_fact = spacetime_lattice::num_pLloops_half(U.getndims() - 1);
    Uij /= double(dims_fact);
    return Uij;
  }

} // namespace operators

// tests done with
//~ double loop = operators::measure_re_arbitrary_loop_lattice(U, {1,1,1,1}, {1,2,1,2},
//{true, true, false, false}); ~ loop += operators::measure_re_arbitrary_loop_lattice(U,
//{1,1,1,1}, {1,3,1,3}, {true, true, false, false}); loop +=
// operators::measure_re_arbitrary_loop_lattice(U, {1,1,1,1}, {2,3,2,3}, {true, true,
// false, false}); std::cout << loop/U.getVolume()*2./U.getndims()/(U.getndims()-1) << "
// "; ~ std::cout << loop/U.getVolume()/2 << " "; ~ loop =
// operators::measure_re_arbitrary_loop_lattice(U, {2,1,2,1}, {0,1,0,1}, {true, true,
// false, false}); loop += operators::measure_re_arbitrary_loop_lattice(U, {1,1,1,1},
// {0,2,0,2}, {true, true, false, false}); loop +=
// operators::measure_re_arbitrary_loop_lattice(U, {1,1,1,1}, {0,3,0,3}, {true, true,
// false, false}); std::cout << loop/U.getVolume()*2./U.getndims()/(U.getndims()-1) <<
// std::endl; ~ std::cout << loop/U.getVolume() << std::endl;
// in measure-u1 agree with results of Wilson-Loops
