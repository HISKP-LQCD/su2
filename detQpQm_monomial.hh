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

// detQpQm monomial_3d : evaluation of det(Q_{+}*Q_{-}) through pseudo-fermions in 2+1 dimensions
template <typename Float, class Group>
class detQpQm_monomial_3d : public monomial<Float, Group> {
public:
Float m0;   // bare mass (in lattice units)

Float TOLERANCE; // tolerance of the CG solver
size_t VERBOSITY; // verbosity of the CG solver
size_t SEED; // seed of the random number generator

// pseudo-fermion field, kept constant along the md trajectory
staggered::spinor_lat_3d<Float, std::complex<Float>> phi;

detQpQm_monomial_3d<Float, Group>(unsigned int _timescale, const Float& m0_val, Float tolerance, const size_t &seed, const size_t& verb)
        : monomial<Float, Group>::monomial(_timescale) {
        m0 = m0_val;

        TOLERANCE = tolerance;
        SEED = seed;
        VERBOSITY = verb;
}


// during the heatbath we generate R and store phi
void heatbath(hamiltonian_field<Float, Group> const &h) override {
        const int N =   h.U->getLt() * h.U->getLx() * h.U->getLy();

        // generating the gaussian R
        const staggered::spinor_lat_3d<Float, std::complex<Float>> R =
                staggered::gaussian_spinor_normalized<Float, std::complex<Float>>(N, SEED);

        // computing the operator Q_{+}
        const staggered::matrix_lat_3d<Float, std::complex<Float>> Qp = staggered::get_Qp_3d<Float, std::complex<Float>>(h.U, (*this).m0);
        (*this).phi = Qp * R;
        monomial<Float, Group>::Hold = R.norm_squared(); // R^{\dagger}*R is real
        return;
}

// after the md trajectory we retain the new R^{\dagger}*R:
// R = Qp^{-1}*phi, where phi is the one computed at the beginning of the trajectory
// and Qp is evaluated from the new gauge config
void accept(hamiltonian_field<Float, Group> const &h) override {
        staggered::matrix_lat_3d<Float, std::complex<Float>> Qp = staggered::get_Qp_3d<Float, std::complex<Float>>(h.U, (*this).m0);// Qp after the update
        staggered::spinor_lat_3d<Float, std::complex<Float>> R = Qp.inv((*this).phi, TOLERANCE, VERBOSITY, SEED); // R = Qp^{-1}*phi
        monomial<Float, Group>::Hnew = R.norm_squared(); // R^{\dagger}*R is real
        return;
}

void derivative(adjointfield<Float, Group> &deriv,
                hamiltonian_field<Float, Group> const &h,
                const Float fac = 1.) const override {
        const int nd = h.U->getndims();
        std::vector<size_t> x(nd), xm(nd), xp(nd); // updated at each loop step
        const size_t Lt = h.U->getLt(), Lx = h.U->getLx(), Ly = h.U->getLy();//, Lz = h.U->getLz();
        const std::vector<size_t> dims = {Lt, Lx, Ly};
        typedef typename accum_type<Group>::type accum;
#pragma omp parallel for
  for(size_t x0 = 0; x0 < Lt; x0++) {
    for(size_t x1 = 0; x1 < Lx; x1++) {
      for(size_t x2 = 0; x2 < Ly; x2++) {
        x = {x0, x1, x2};
        xm = x; xp = x;
        for(size_t mu = 0; mu < h.U->getndims(); mu++) {
        xm[mu] -= 1; // x - mu
        xp[mu] += 1; // x + mu
        
        accum S = 
        conj(phi(x, dims)) * phi(xm, dims) 
        + get_deriv<double>(h.U(x, mu).dagger()) * conj(phi(xp, dims)) * phi(x, dims);
        deriv(x, mu) += (1.0/2.0)*staggered::eta(x, mu)*S;
        
        xm[mu] += 1; // =x again
        xp[mu] -= 1; // =x again
        }
      }
    }
  }
        return;
}

};
