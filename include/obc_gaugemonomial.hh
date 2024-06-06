/**
 * @file obc-gaugemonomial.hh
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief gauge monomial S_G with open boundary conditions
 * @version 0.1
 * @date 2022-01-30
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once
#include "adjointfield.hh"
#include "gauge_energy.hpp"
#include "gaugeconfig.hh"
#include "geometry.hh"
#include "get_staples.hh"
#include "hamiltonian_field.hh"
#include "monomial.hh"
#include "obc_weights.hh"
#include "su3.hh"
#include "su3_accum.hh"
#include "su2.hh"
#include "u1.hh"

#include <complex>
#include <vector>
#include <cmath>

namespace obc { // open boundary conditions

/**

 * @brief sum of the real trace of all gauge links with periodic boundary conditions 
 * 
 * @tparam T gauge group
 * @param U gauge configuration
 * @return double 
 */
  template <class T>
  double retr_sum_realtrace(const gaugeconfig<T> &U){
    double res = 0.;
  #pragma omp parallel for reduction(+ : res)
    for (size_t x0 = 0; x0 < U.getLt(); x0++){
      for (size_t x1 = 0; x1 < U.getLx(); x1++){
        for (size_t x2 = 0; x2 < U.getLy(); x2++){
          for (size_t x3 = 0; x3 < U.getLz(); x3++){
            const std::vector<size_t> x = {x0, x1, x2, x3};
            for(size_t mu = 0; mu < U.getndims() -1; mu++){
              res += retrace(U(x, mu));
            }
          }
        }
      }
    }
  return res;
  }
  /**
   * @brief real part of the trace of the Wilson plaquette (with spatial open boundary
   * conditions)
   *
   *
   * @tparam T gauge group
   * @param U gauge configuration
   * @param spatial_only: true when only the plaquettes with mu, nu > 0 are calculated
   * @param xi anisotropy in the action: S_G \supset (\beta/xi)*P_{0i} + (xi*\beta)*P_{ij}
   * @return double
   */
  template <class T>
  double retr_sum_Wplaquettes(const gaugeconfig<T> &U,
                              const obc::weights &w,
                              const double &xi = 1.0,
                              const bool &anisotropic = false,
                              const bool &spatial_only = false) {
    double res = 0.;
    const size_t startmu = spatial_only; // 0 if spatial_only==false, 1 if spatial_only==true

    if (!anisotropic) {
#pragma omp parallel for reduction(+ : res)
      for (size_t x0 = 0; x0 < U.getLt(); x0++) {
        for (size_t x1 = 0; x1 < U.getLx(); x1++) {
          for (size_t x2 = 0; x2 < U.getLy(); x2++) {
            for (size_t x3 = 0; x3 < U.getLz(); x3++) {
              const std::vector<size_t> x = {x0, x1, x2, x3};
              const double wx = w(x);
              std::vector<size_t> xplusmu = x;
              std::vector<size_t> xplusnu = x;
              for (size_t mu = startmu; mu < U.getndims() - 1; mu++) {
                for (size_t nu = mu + 1; nu < U.getndims(); nu++) {
                  xplusmu[mu] += 1;
                  xplusnu[nu] += 1;

                  const double wxplusmu = w(xplusmu);
                  const double wxplusnu = w(xplusnu);
                  res += wx * wxplusmu * wxplusnu *
                         retrace(U(x, mu) * U(xplusmu, nu) * U(xplusnu, mu).dagger() *
                                 U(x, nu).dagger());

                  xplusmu[mu] -= 1;
                  xplusnu[nu] -= 1;
                }
              }
            }
          }
        }
      }
    }
    // 2n option - anisotropic lattice present
    if (anisotropic) {
#pragma omp parallel for reduction(+ : res)
      for (size_t x0 = 0; x0 < U.getLt(); x0++) {
        for (size_t x1 = 0; x1 < U.getLx(); x1++) {
          for (size_t x2 = 0; x2 < U.getLy(); x2++) {
            for (size_t x3 = 0; x3 < U.getLz(); x3++) {
              const std::vector<size_t> x = {x0, x1, x2, x3};
              const double wx = w(x);
              std::vector<size_t> xplusmu = x;
              std::vector<size_t> xplusnu = x;
              for (size_t mu = startmu; mu < U.getndims() - 1; mu++) {
                for (size_t nu = mu + 1; nu < U.getndims(); nu++) {
                  xplusmu[mu] += 1;
                  xplusnu[nu] += 1;

                  double eta = xi;
                  if ((mu == 0) ^ (nu == 0)) {
                    // at least one direction is temporal (but not both)
                    eta = 1.0 / xi;
                  }
                  const double wxplusmu = w(xplusmu);
                  const double wxplusnu = w(xplusnu);
                  res += eta * wx * wxplusmu * wxplusnu *
                         retrace(U(x, mu) * U(xplusmu, nu) * U(xplusnu, mu).dagger() *
                                 U(x, nu).dagger());

                  xplusmu[mu] -= 1;
                  xplusnu[nu] -= 1;
                }
              }
            }
          }
        }
      }
    }
    return res;
  }

  /**
   * @brief Get the S_G action for the HMC
   * Get the action subtracted from unnecessary constant factors (unphisical).
   *
   * Returns:
   * (\beta/Nc)*\sum_{x} \sum_{\mu<\nu} \eta_{\mu \nu} Re(Tr(P_{\mu \nu}(x)))
   * with \eta_{0i}=1/\xi and \eta_{ij}=\xi
   *
   * where beta = 2*N_c/g_0^2
   *
   * @tparam Float
   * @tparam Group
   * @param U gauge configuration
   * @param xi anisotropy (==1.0 when the lattice is isotropic)
   * @param anisotropic true when considering the anisotropy
   * @return Float
   */
  template <typename Float, class Group>
  Float get_S_G_hmc(const gaugeconfig<Group> &U,
                    const obc::weights &w,
                    const double &xi = 1.0,
                    const bool &anisotropic = false) {
    const double fact = U.getBeta() / double(U.getNc());
    const Float res = -fact * obc::retr_sum_Wplaquettes(U, w, xi, anisotropic);
    return res;
  }

  /**
   * @brief Full gauge action
   * Returns
   *
   * S_G = (\beta/Nc)*\sum_{x} \sum_{\mu<\nu}
   * \eta_{\mu \nu} Re(Tr(1 - P_{\mu\nu}(x)))
   *
   * where beta = 2*N_c/g_0^2
   *
   * @tparam Float
   * @tparam Group
   * @param U
   * @param xi
   * @param anisotropic
   * @return Float
   */
  template <typename Float, class Group>
  Float get_S_G(const gaugeconfig<Group> &U,
                const obc::weights &w,
                const double &xi = 1.0,
                const bool &anisotropic = false) {
    const size_t Dims_fact = spacetime_lattice::num_pLloops_half(U.getndims());
    const Float res = U.getBeta() * (U.getVolume() * Dims_fact) +
                      obc::get_S_G_hmc<Group>(U, xi, anisotropic);
    return res;
  }

  // gauge monomial
  template <typename Float, class Group>
  class gaugemonomial : public monomial<Float, Group> {
  public:
    gaugemonomial<Float, Group>(unsigned int _timescale,
                                const obc::weights &_w,
                                const double &_xi = 1.0)
      : monomial<Float, Group>::monomial(_timescale) {
      w = _w;
      xi = _xi;
      if (_xi != 1.0) {
        anisotropic = true;
      }
    }

    // S_g = sum_x sum_{mu<nu} beta*(1- 1/Nc*Re[Tr[U_{mu nu}]])
    // beta = 2*N_c/g_0^2
    void heatbath(hamiltonian_field<Float, Group> const &h) override {
      monomial<Float, Group>::Hold =
        obc::get_S_G_hmc<Float, Group>(*h.U, (*this).w, (*this).xi, (*this).anisotropic);
      return;
    }
    void accept(hamiltonian_field<Float, Group> const &h) override {
      monomial<Float, Group>::Hnew =
        obc::get_S_G_hmc<Float, Group>(*h.U, (*this).w, (*this).xi, (*this).anisotropic);
      return;
    }

    /**
     * @brief updating the derivative of the action with respect to the gauge field
     * contribution in this monomial
     *
     * @param deriv derivative to be updated
     * @param h hamiltonian field (Hybrid Monte Carlo evolution)
     * @param fac multiplicative facator of the derivative (1.0 for the HMC but != 1 in
     * other cases e.g. the gradient flow evolution)
     */
    void derivative(adjointfield<Float, Group> &deriv,
                    hamiltonian_field<Float, Group> const &h,
                    const Float fac = 1.) const override {
      typedef typename accum_type<Group>::type accum;
      const double num_fact_i = fac * h.U->getBeta() / double(h.U->getNc());
#pragma omp parallel for
      for (size_t x0 = 0; x0 < h.U->getLt(); x0++) {
        for (size_t x1 = 0; x1 < h.U->getLx(); x1++) {
          for (size_t x2 = 0; x2 < h.U->getLy(); x2++) {
            for (size_t x3 = 0; x3 < h.U->getLz(); x3++) {
              const std::vector<size_t> x = {x0, x1, x2, x3};
              const double w_x = (*this).w(x);
              std::vector<size_t> xpmu = x, xpnu = x, xmnu = x;
              for (size_t mu = 0; mu < h.U->getndims(); mu++) {
                xpmu[mu]++; // x + mu

                const double w_xpmu = (*this).w(xpmu);
                accum S;
                for (size_t nu = 0; nu < h.U->getndims(); nu++) {
                  if (mu == nu) {
                    continue;
                  }
                  xpnu[nu]++; // x + nu
                  xmnu[nu]--; // x - nu

                  const double w_xpnu = (*this).w(xpnu);
                  const accum S1 = w_x * w_xpmu * w_xpnu *
                                   get_staple_up<accum>(*h.U, mu, nu, xpmu, xpnu, x);
                  S += S1;

                  xpmu[nu]--; // x + mu - nu
                  const double w_xmnu = (*this).w(xmnu);
                  const double w_xpmumnu =
                    (*this).w(xpmu); // Here: xpmu == x + mu - nu !!!
                  const accum S2 = w_x * w_xmnu * w_xpmumnu *
                                   get_staple_down<accum>(*h.U, mu, nu, xpmu, xmnu);
                  S += S2;
                  xpmu[nu]++; // x + mu

                  xpnu[nu]--; // x
                  xmnu[nu]++; // x
                }
                S = (*h.U)(x, mu) * S; // U*A in eq. 8.40 in Gattringer&Lang
                deriv(x, mu) += num_fact_i * get_deriv(S);

                xpmu[mu]--; // x
              }
            }
          }
        }
      }
      return;
    }

  private:
    // size_t Dims_fact; // d*(d-1)/2

    bool anisotropic = false;
    double xi; // bare anisotropy xi
    obc::weights w; // weights for obc
  };

} // namespace obc
