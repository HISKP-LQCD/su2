// detDdagD_monomial
// Evaluation of det(D*D^{\dagger}) with the pseudo-fermion action 
// as in eq. (10) of https://www.sciencedirect.com/science/article/pii/0550321389903246
// or s eq. (8.9) of https://link.springer.com/book/10.1007/978-3-642-01850-3

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

// detDdagD monomial_4d : evaluation of det(D^{\dagger}*D) through pseudo-fermions
template <typename Float, class Group>
class detDdagD_monomial_4d : public monomial<Float, Group> {
typedef std::complex<Float> Complex;

public:
Float m0;   // bare mass (in lattice units)

Float TOLERANCE; // tolerance of the CG solver
size_t VERBOSITY; // verbosity of the CG solver
size_t SEED; // seed of the random number generator

// pseudo-fermion field phi, kept constant along the MD trajectory
staggered::spinor_lat_4d<Float, Complex> phi;

detDdagD_monomial_4d<Float, Group>(unsigned int _timescale, const Float& m0_val, Float tolerance, const size_t &seed, const size_t& verb)
        : monomial<Float, Group>::monomial(_timescale) {
        m0 = m0_val;

        TOLERANCE = tolerance;
        SEED = seed;
        VERBOSITY = verb;
}


// during the heatbath we generate R and store phi
void heatbath(hamiltonian_field<Float, Group> const &h) override {
  const int N = h.U->getLt() * h.U->getLx() * h.U->getLy() * h.U->getLz(); // total number of lattice points
  
  // generating the gaussian R
  const staggered::spinor_lat_4d<Float, Complex> R =
  staggered::gaussian_spinor_normalized<Float, Complex>(N, 0.0, 1.0 / sqrt(2), SEED); // e^{-x^2} has sigma=1/sqrt(2)
  // applying the operator D^{\dagger} to R
  (*this).phi = staggered::apply_Ddag<Float, Complex, Group>(h.U, (*this).m0, R);
  monomial<Float, Group>::Hold = R.norm_squared(); // R^{\dagger}*R is real
  return;
}


/**
 * After the MD trajectory we compute the new R^{\dagger}*R : 
 * R = D*chi = D * ( (D^{\dagger}*D)^{-1} * phi ), 
 * where phi is the one computed at the beginning of the trajectory
 * and D, D^{\dagger} are evaluated from the new gauge config  
 */
void accept(hamiltonian_field<Float, Group> const &h) override {
  // Operator D^{\dagger}*D . Hermitean and invertible -> can apply the CG inversion

  // const staggered::matrix_lat_4d<Float, Complex> DdagD = 
  //   staggered::get_DdagD_4d<Float, Complex>(h.U, (*this).m0);
  const staggered::DdagD_matrix_lat<Float, Complex, Group> DdagD(h.U, (*this).m0);

  // applying Ddag*D to \chi
  const staggered::spinor_lat_4d<Float, Complex> chi = 
    DdagD.inv((*this).phi, TOLERANCE, VERBOSITY, SEED);

  const staggered::spinor_lat_4d<Float, Complex> R = 
    staggered::apply_D(h.U, (*this).m0, chi);
  monomial<Float, Group>::Hnew = R.norm_squared(); // R^{\dagger}*R is real
  return;
}


/**
 * S_F  = \phi^{\dagger} * M^{-1} * \phi 
 * dS_F = \phi^{\dagger} * M^{-1} * dM * M^{-1} * \phi = \chi^{\dagger} * dM * \chi  
 * where \chi = M^{-1} * \phi
 * see eq. (8.44) of Gattringer&Lang
 */
void derivative(adjointfield<Float, Group> &deriv,
                hamiltonian_field<Float, Group> const &h,
                const Float fac = 1.) const override {
  typedef typename accum_type<Group>::type accum;

  const staggered::DdagD_matrix_lat<Float, Complex, Group> DdagD(h.U, (*this).m0);
  const staggered::spinor_lat_4d<Float, Complex> chi = 
    DdagD.inv((*this).phi, TOLERANCE, VERBOSITY, SEED);

  const size_t Lt = h.U->getLt(), Lx = h.U->getLx(), Ly = h.U->getLy(), Lz = h.U->getLz();
  const std::vector<size_t> dims = {Lt, Lx, Ly, Lz}; // vector of dimensions
  const size_t nd = h.U->getndims();

#pragma omp parallel for
  for(size_t x0 = 0; x0 < Lt; x0++) {
    for(size_t x1 = 0; x1 < Lx; x1++) {
      for(size_t x2 = 0; x2 < Ly; x2++) {
        for(size_t x3 = 0; x3 < Lz; x3++) {
          const std::vector<size_t> x = {x0, x1, x2, x3};
          // std::vector<size_t> xm = x; // will be x+mu
          // std::vector<size_t> xp = x; // sill be x-mu
          // std::vector<size_t> xpm = x; // will be x + mu - nu
          // std::vector<size_t> xmm = x; // will be x - mu - nu
          // std::vector<size_t> xmp = x; // will be  x - mu + nu
          // std::vector<size_t> xpp = x; // will be  x + mu + nu

          // accum conj_chi_x = conj(chi(x, dims)); 
          
          accum D;

          for(size_t mu = 0; mu < nd; mu++) {
          //   xp[mu] += 1; // x + mu
          //   xm[mu] -= 1; // x - mu
          //   xpp[mu] += 1; // see later in the loop
          //   xpm[mu] += 1; // see later in the loop
          //   xmp[mu] -= 1; // see later in the loop
          //   xmm[mu] -= 1; // see later in the loop
          //   //
          //   const Float fact_mu   = (1.0/2.0) * staggered::eta(x, mu);
          //   for(size_t nu = 0; nu < nd; nu++) {
          //     xpp[mu] += 1; // x+mu+nu
          //     xpm[mu] -= 1; // x+mu-nu
          //     xmp[mu] += 1; // x-mu+nu
          //     xmm[mu] -= 1; // x-mu-nu

          //     const Float fact_nu_p = (1.0/2.0) * staggered::eta(xp, nu);
          //     const Float fact_nu_m = (1.0/2.0) * staggered::eta(xm, nu);

          //     const Float fact1 = fact_mu * fact_nu_p;
          //     const Float fact2 = fact_mu * fact_nu_m;

          //     accum D11, D12, D21, D22; 

          //     D11 = + conj_chi_x * fact1 * (*h.U)(x, mu).dagger() * (*h.U)(xp, mu) * chi(xpp, dims);
          //     D += D11;
              
          //     D12 = - conj_chi_x * fact1 * (*h.U)(x, mu).dagger() * (*h.U)(xpm, mu).dagger() * chi(xmm, dims);
          //     D += D12;

          //     D21 = - conj_chi_x * fact2 * (*h.U)(xm, mu) * (*h.U)(xm, mu) * chi(xmm, dims);              
          //     D += D21;

          //     D22 = + conj_chi_x * fact2 * (*h.U)(xm, mu) * (*h.U)(xmm, mu).dagger() * chi(xmm, dims);
          //     D += D22;

          //     xpp[mu] -= 1; // x+mu again
          //     xpm[mu] += 1; // x+mu again
          //     xmp[mu] -= 1; // x-mu again
          //     xmm[mu] += 1; // x-mu again
          //   }

          //   accum E1 = + conj_chi_x * (*this).m0 * fact_mu * ( (*h.U)(x, mu) + (*h.U)(x, mu).dagger() ) * chi(xp, dims);
          //   D += E1;

          //   accum E2 = - conj_chi_x * (*this).m0 * fact_mu * ( (*h.U)(xm, mu).dagger() + (*h.U)(xm, mu) ) * chi(xm, dims);
          //   D += E2;

          const staggered::spinor_lat_4d<Float, Complex> chi1 = staggered::apply_Ddag(h.U, (*this).m0, chi);
          const staggered::spinor_lat_4d<Float, Complex> chi2 = staggered::apply_der_Ddag(x, mu, h.U, (*this).m0, chi);
          const staggered::spinor_lat_4d<Float, Complex> chi_prime = 
            staggered::apply_der_D(x, mu, h.U, (*this).m0, chi1) +
            staggered::apply_D(h.U, (*this).m0, chi2);

          // see eq. (12) of https://www.sciencedirect.com/science/article/pii/0550321389903246
          D = staggered::complex_dot_product(chi, chi_prime);
          D = 2.0 * (*h.U)(x, mu) * D;
          deriv(x, mu) += get_deriv<Float>(D); 

            // xm[mu] += 1; // =x again
            // xp[mu] -= 1; // =x again
            // xpp[mu] -= 1; // =x again
            // xpm[mu] -= 1; // =x again
            // xmp[mu] += 1; // =x again
            // xmm[mu] += 1; // =x again
          }
        }
      }
    }
  }
  return;
}

};
