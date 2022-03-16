// detDDdag_monomial
/*
  Simone Romiti - simone.romiti@uni-bonn.de

  Evaluation of det(D*D^{\dagger}) with the pseudo-fermion action
  as in eq. (10) of https://www.sciencedirect.com/science/article/pii/0550321389903246
  or s eq. (8.9) of https://link.springer.com/book/10.1007/978-3-642-01850-3
*/

#pragma once

#include "adjointfield.hh"
#include "gauge_energy.hh"
#include "gaugeconfig.hh"
#include "get_staples.hh"
#include "hamiltonian_field.hh"
#include "monomial.hh"
#include "su2.hh"
#include "u1.hh"
#include <complex>
#include <random>
#include <vector>

#include "staggered.hpp" // spinor object

// detDDdag monomial : evaluation of det(D*D^{\dagger}) through pseudo-fermions
template <typename Float, class Group>
class detDDdag_monomial : public monomial<Float, Group> {
  typedef std::complex<Float> Complex;

public:
  Float m0; // bare mass (in lattice units)

  Float TOLERANCE; // tolerance of the CG solver
  size_t VERBOSITY; // verbosity of the CG solver
  size_t SEED; // seed of the random number generator

  // pseudo-fermion field phi, kept constant along the MD trajectory
  staggered::spinor_lat<Float, Complex> phi;

  detDDdag_monomial<Float, Group>(unsigned int _timescale,
                                     const Float &m0_val,
                                     const Float &tolerance,
                                     const size_t &seed,
                                     const size_t &verb)
    : monomial<Float, Group>::monomial(_timescale) {
    m0 = m0_val;

    TOLERANCE = tolerance;
    SEED = seed;
    VERBOSITY = verb;
  }

  // during the heatbath we generate R and store phi
  void heatbath(hamiltonian_field<Float, Group> const &h) override {
    const int N = h.U->getVolume(); // total number of lattice points

    const size_t Lt = h.U->getLt(), Lx = h.U->getLx(), Ly = h.U->getLy(),
                 Lz = h.U->getLz();
    const std::vector<size_t> dims = {Lt, Lx, Ly, Lz}; // vactor of spacetime dimensions

    // generating the gaussian R
    const staggered::spinor_lat<Float, Complex> R =
      staggered::gaussian_spinor_normalized<Float, Complex>(
        dims, N, 0.0, 1.0 / sqrt(2), SEED); // e^{-x^2} has sigma=1/sqrt(2)
    // applying the operator D^{\dagger} to R
    (*this).phi = staggered::apply_D<Float, Complex, Group>(h.U, (*this).m0, R);

    monomial<Float, Group>::Hold = R.norm_squared(); // R^{\dagger}*R is real
    return;
  }

  /**
   * After the MD trajectory we compute the new R^{\dagger}*R :
   * R = D*chi = D * ( (D*D^{\dagger})^{-1} * phi ),
   * where phi is the one computed at the beginning of the trajectory
   * and D, D^{\dagger} are evaluated from the new gauge config
   */
  void accept(hamiltonian_field<Float, Group> const &h) override {
    // Operator D*D^{\dagger} . Hermitian and invertible -> can apply the CG inversion

    const staggered::DDdag_matrix_lat<Float, Complex, Group> DDdag(h.U, (*this).m0);

    // applying Ddag*D to \chi
    const staggered::spinor_lat<Float, Complex> chi =
      DDdag.inv((*this).phi, TOLERANCE, VERBOSITY, SEED);

    const staggered::spinor_lat<Float, Complex> R =
      staggered::apply_Ddag(h.U, (*this).m0, chi);

    monomial<Float, Group>::Hnew = R.norm_squared(); // R^{\dagger}*R is real
    return;
  }

  /**
   * S_F  = \phi^{\dagger} * M^{-1} * \phi
   * dS_F = \phi^{\dagger} * M^{-1} * dM * M^{-1} * \phi = \chi^{\dagger} * dM * \chi
   * where \chi = M^{-1} * \phi
   * and M = D*D^{\dagger}
   * see eq. (8.44) of Gattringer&Lang
   */
  void derivative(adjointfield<Float, Group> &deriv,
                  hamiltonian_field<Float, Group> const &h,
                  const Float fac = 1.) const override {
    typedef typename accum_type<Group>::type accum;

    const staggered::DDdag_matrix_lat<Float, Complex, Group> DDdag(h.U, (*this).m0);

    const staggered::spinor_lat<Float, Complex> chi =
      DDdag.inv((*this).phi, TOLERANCE, VERBOSITY, SEED);


    const size_t Lt = h.U->getLt(), Lx = h.U->getLx(), Ly = h.U->getLy(),
                 Lz = h.U->getLz();
    const std::vector<size_t> dims = {Lt, Lx, Ly, Lz}; // vector of dimensions
    const size_t nd = h.U->getndims();

    const staggered::spinor_lat<Float, Complex> chi1 =
      staggered::apply_Ddag(h.U, (*this).m0, chi);

    const Complex i(0.0, 1.0);

#pragma omp parallel for
    for (int x0 = 0; x0 < Lt; x0++) {
      for (int x1 = 0; x1 < Lx; x1++) {
        for (int x2 = 0; x2 < Ly; x2++) {
          for (int x3 = 0; x3 < Lz; x3++) {
            const std::vector<int> x = {x0, x1, x2, x3};
            std::vector<int> xm = x, xp = x;
            for (size_t mu = 0; mu < nd; mu++) {
              xm[mu]--; // x - mu
              xp[mu]++; // x + mu

              const Float eta_x_mu = staggered::eta(x, mu);
              const Float eta_xp_mu = staggered::eta(xp, mu);

              const Complex v_xp =
                (1.0 / 2.0) * eta_x_mu * (+i) * (*h.U)(x, mu) * chi1(xp);
              const Complex v_x =
                -(1.0 / 2.0) * eta_xp_mu * (-i) * (*h.U)(x, mu).dagger() * chi1(x);

              // derivative of S_F with respect to U_{\mu}(x)
              accum derSF = conj(chi(x)) * v_xp + conj(chi(xp)) * v_x;
              derSF *= -1.0;

              derSF *= i; // get_deriv gives the imaginary part: Re(z) = Im(i*z)
              deriv(x, mu) += get_deriv<Float>(derSF);

              xm[mu]++; // =x again
              xp[mu]--; // =x again
            }
          }
        }
      }
    }
    return;
  }
};
