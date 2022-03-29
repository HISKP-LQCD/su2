// staggered.hpp
/**
 Simone Romiti - simone.romiti@uni-bonn.de

 Staggered fermions routines
*/

#pragma once

#include <array>
#include <complex>
#include <random>
#include <vector>

#include "adjointfield.hh"
#include "gauge_energy.hh"
#include "gaugeconfig.hh"
#include "get_staples.hh"
#include "hamiltonian_field.hh"
#include "include/geometry.hh"
#include "monomial.hh"
#include "su2.hh"
#include "u1.hh"

#include "solver_type.hh" // generic type of solver

// given A(matrix) and b(vector), find x = A^{-1}*b
#include "CG.hpp" // standard conjugate gradient 
#include "BiCGStab.hpp" // Bi Conjugate Gradient

// std::vector<double> g_vec_eta_x_mu; // values of \eta_{\mu}(x)

namespace staggered {

  // eta_{\mu}(x) as in eq. (16) of https://arxiv.org/pdf/2112.14640.pdf
  // or eq. (8) of https://www.sciencedirect.com/science/article/pii/0550321389903246
  double eta(const std::vector<int> &x, const size_t &mu) {
    double s = 0;
// TO BE TESTED #pragma omp parallel for reduction(+: s)
    for (int nu = 0; nu < mu; nu++) {
      s += x[nu];
    }
    return std::pow(-1.0, s);
  }

//   /* fills the array arr_eta */
//   void generate_eta_x_mu(const std::vector<size_t>& dims, const size_t &ndims) {
//     const size_t Lt = dims[0], Lx = dims[1], Ly = dims[2], Lz = dims[3];
//     size_t N = Lt*Lx*Ly*Lz*ndims;
//     g_vec_eta_x_mu.resize(N);
//     const geometry_d g(Lx, Ly, Lz, Lt, ndims); // note the order
// #pragma omp parallel for
//     for (int x0 = 0; x0 < Lt; x0++) {
//       for (int x1 = 0; x1 < Lx; x1++) {
//         for (int x2 = 0; x2 < Ly; x2++) {
//           for (int x3 = 0; x3 < Lz; x3++) {
//             const std::vector<int> x = {x0, x1, x2, x3};
//             for (size_t mu = 0; mu < ndims; mu++) {
//               const size_t i = g.getIndex(x0, x1, x2, x3);
//               g_vec_eta_x_mu[i] = eta(x, mu);
//             }
//           }
//         }
//       }
//     }
//     return;
//   }

  // index of the lattice point given the dimensions
  size_t txyz_to_index(const size_t &t,
                       const size_t &x,
                       const size_t &y,
                       const size_t &z,
                       const size_t &Lt,
                       const size_t &Lx,
                       const size_t &Ly,
                       const size_t &Lz) {
    const geometry g(Lx, Ly, Lz, Lt); // note the order
    return g.getIndex(t, x, y, z);
  }

  // overload : x={x0,x1,x2,x3}, dims={Lt,Lx,Ly,Lz}
  size_t txyz_to_index(const std::vector<int> &x, const std::vector<size_t> &dims) {
    return txyz_to_index(x[0], x[1], x[2], x[3], dims[0], dims[1], dims[2], dims[3]);
  }

  // vector of staggered "spinors" (no Dirac structure) for all the points of the lattice
  template <class Float, class Type> class spinor_lat {
  private:
    std::vector<Type> Psi;
    std::vector<size_t> dims; // spacetime dimensions : {Lt, Lx, Ly, Lz}

  public:
    spinor_lat() {}
    ~spinor_lat() {}

    Type &operator()(const size_t &i) { return Psi[i]; }
    Type operator()(const size_t &i) const { return Psi[i]; }

    Type &operator()(const std::vector<int> &x) {
      const size_t i = txyz_to_index(x, dims);
      return (*this)[i];
    }

    Type operator()(const std::vector<int> &x) const {
      const size_t i = txyz_to_index(x, dims);
      return (*this)[i];
    }

    Type &operator[](const size_t &i) { return Psi[i]; }
    Type operator[](const size_t &i) const { return Psi[i]; }

    spinor_lat(const std::vector<size_t> &_dims, const size_t &n) {
      Psi.resize(n);
      dims = _dims;
    }

    spinor_lat(const std::vector<size_t> &_dims, const size_t &n, const Type &val) {
      std::vector<Type> v(n, val);
      (*this).Psi = v;
      dims = _dims;
    }

    std::vector<size_t> get_dims() const { return dims; }

    size_t size() const { return Psi.size(); }
    void resize(const size_t &n, const Type &val = (Type)0.0) { Psi.resize(n, val); }

    spinor_lat<Float, Type> operator/(const Type &lambda) {
      const int N = (*this).size();
      spinor_lat<Float, Type> phi(this->get_dims(), N);
      for (size_t i = 0; i < N; i++) {
        phi[i] = Psi[i] / lambda;
      }
      return phi;
    }

    void operator+=(const spinor_lat<Float, Type> &psi) { (*this) = (*this) + psi; }

    Type dot(const spinor_lat& b) const{
      return complex_dot_product((*this), b); // (*this)^{\dagger} * b
    }

    Float norm_squared() const {
      return complex_dot_product((*this), (*this)).real(); // A^{\dagger} * A is real
    }

    Float norm() const { return sqrt(this->norm_squared()); }

  }; // class spinor_lat

  // a + lambda*b
  template <class Float, class Type, class Type_lambda>
  spinor_lat<Float, Type> a_plus_lambda_b(const spinor_lat<Float, Type> &a,
                                          const Type_lambda &lambda,
                                          const spinor_lat<Float, Type> &b) {
    const int N = a.size();
    spinor_lat<Float, Type> c(a.get_dims(), N);
    for (size_t i = 0; i < N; i++) {
      c[i] = a[i] + lambda * b[i];
    }
    return c;
  }

  template <class Float, class Type>
  spinor_lat<Float, Type> operator+(const spinor_lat<Float, Type> &a,
                                    const spinor_lat<Float, Type> &b) {
    return a_plus_lambda_b(a, 1.0, b);
  }

  template <class Float, class Type>
  spinor_lat<Float, Type> operator-(const spinor_lat<Float, Type> &a,
                                    const spinor_lat<Float, Type> &b) {
    return a_plus_lambda_b(a, -1.0, b);
  }

  // change of sign : psi --> -psi
  template <class Float, class Type>
  spinor_lat<Float, Type> operator-(const spinor_lat<Float, Type> &psi) {
    const spinor_lat<Float, Type> v(psi.get_dims(), psi.size());
    return (v - psi);
  }

  template <class Float, class Type>
  spinor_lat<Float, Type> gaussian_spinor_normalized(const std::vector<size_t> &dims,
                                                     const size_t &n,
                                                     const Float &avr,
                                                     const Float &sigma,
                                                     const size_t &seed) {
    std::normal_distribution<Float> dis{avr, sigma};
    std::mt19937 gen_re(seed), gen_im(seed + 1);
    spinor_lat<Float, Type> psi_gauss(dims,
                                      n); // spacetime dimensions are irrelevant here
    Type norm2 = 0.0;
    for (size_t i = 0; i < n; i++) { // lattice points
      Type x = dis(gen_re);
      psi_gauss[i] = x;
      norm2 += (conj(x) * x).real(); // x^{\dagger}*x is real
    }
    return psi_gauss / sqrt(norm2);
  }

  // \sum_{i} A_i^{\dagger}*B_i
  template <class Float, class Type>
  Type complex_dot_product(const spinor_lat<Float, Type> &A,
                           const spinor_lat<Float, Type> &B) {
    const int N = A.size();
    Type sum = 0.0;

    // #pragma omp parallel for
    for (size_t i = 0; i < N; i++) {
      sum += conj(A[i]) * B[i];
    }

    return sum;
  }

  template <class Float, class Type>
  Type operator*(const spinor_lat<Float, Type> &A, const spinor_lat<Float, Type> &B) {
    return complex_dot_product(A, B);
  }

  template <class Float, class Type>
  spinor_lat<Float, Type> operator*(const Type &lambda,
                                    const spinor_lat<Float, Type> &psi) {
    const spinor_lat<Float, Type> v(psi.get_dims(), psi.size());
    return a_plus_lambda_b(v, lambda, psi);
  }

  template <class Float, class Type> Float norm(const spinor_lat<Float, Type> &psi) {
    return psi.norm();
  }

  // class for D^{\dagger}*D
  template <class Float, class Type, class Group> class DDdag_matrix_lat {
  public:
    gaugeconfig<Group> *U;
    Float m;

    DDdag_matrix_lat() {}
    ~DDdag_matrix_lat() {}

    DDdag_matrix_lat(gaugeconfig<Group> *_U, const Float &_m) {
      U = _U;
      m = _m;
    }

    // rows = cols = N = number of lattice points
    size_t rows() const { return U->getVolume(); }
    size_t cols() const { return this->rows(); }

    spinor_lat<Float, Type> inv(const spinor_lat<Float, Type> &psi,
                                const std::string& solver,
                                const Float &tol,
                                const size_t &verb,
                                const size_t &seed) const {
      typedef spinor_lat<Float, Type> LAvector;

      typedef DDdag_matrix_lat<Float, Type, Group> LAmatrix;

      typedef CG::LinearCG<Float, Type, LAmatrix, LAvector> cg;        
      typedef typename solver_type<cg>::type svr_type;
      if (solver=="BiCGStab"){
        typedef BiCGStab::LinearBiCGStab<Float, Type, LAmatrix, LAvector> bcgstab;
        typedef typename solver_type<bcgstab>::type svr_type;
      }

      svr_type SVR((*this), psi);

      const size_t N = psi.size();

      const LAvector phi0 = staggered::gaussian_spinor_normalized<Float, Complex>(
        psi.get_dims(), N, 0.0, 10.0, seed);

      if (verb > 1) {
        std::cout << "Calling the CG solver.\n";
      }

      SVR.solve(phi0, tol, verb);
      return SVR.get_solution();
    }
  };

  /**
   * @brief (D*D^{\dagger})*psi
   *
   * @tparam Float : double
   * @tparam Type : std::complex<double>
   * @tparam Group : su1 or su2
   * @param M : wrapper for D*D^{\dagger} matrix.
   *  Contains the information only about the configuration U and the quark mass
   * @param psi : spinor to which we apply D*D^{\dagger}
   * @return spinor_lat<Float, Type>
   */
  template <class Float, class Type, class Group>
  spinor_lat<Float, Type> operator*(const DDdag_matrix_lat<Float, Type, Group> &M,
                                    const spinor_lat<Float, Type> &psi) {
    gaugeconfig<Group> *U = M.U;
    const Float m = M.m;

    const size_t Lt = U->getLt(), Lx = U->getLx(), Ly = U->getLy(), Lz = U->getLz();
    const size_t nd = U->getndims();
    const std::vector<size_t> dims = psi.get_dims(); // vector of dimensions

    const int N = psi.size();
    spinor_lat<Float, Type> phi = apply_D(U, m, apply_Ddag(U, m, psi));
//  spinor_lat<Float, Type> phi(dims, N);
// #pragma omp parallel for
//     for (int x0 = 0; x0 < Lt; x0++) {
//       for (int x1 = 0; x1 < Lx; x1++) {
//         for (int x2 = 0; x2 < Ly; x2++) {
//           for (int x3 = 0; x3 < Lz; x3++) {
//             const std::vector<int> x = {x0, x1, x2, x3};
//             std::vector<int> xm = x; // will be x-mu
//             std::vector<int> xp = x; // will be x+mu
//             std::vector<int> xpp = x; // will be  x + mu + nu
//             std::vector<int> xpm = x; // will be x + mu - nu
//             std::vector<int> xmp = x; // will be  x - mu + nu
//             std::vector<int> xmm = x; // will be x - mu - nu
//             for (size_t mu = 0; mu < nd; mu++) {
//               xp[mu]++; // x + mu
//               xm[mu]--; // x - mu

//               xpp[mu]++; // see later in the loop
//               xpm[mu]++; // see later in the loop
//               xmp[mu]--; // see later in the loop
//               xmm[mu]--; // see later in the loop

//               const Float fact_eta_x_mu = (1.0 / 2.0) * eta(x, mu);
//               const Float fact_eta_xp_mu = (1.0 / 2.0) * eta(xp, mu);
//               const Float fact_eta_xm_mu = (1.0 / 2.0) * eta(xm, mu);

//               for (size_t nu = 0; nu < nd; nu++) {
//                 xpp[nu]++; // x+mu+nu
//                 xpm[nu]--; // x+mu-nu
//                 xmp[nu]++; // x-mu+nu
//                 xmm[nu]--; // x-mu-nu

//                 // std::cout << "mu " << mu << " nu " << nu << "\n";
//                 // std::cout << "x: " << x[0] << "," << x[1] << "," << x[2] << "," << x[3]
//                 //           << "\n";
//                 // std::cout << "xpp: " << xpp[0] << "," << xpp[1] << "," << xpp[2] << ","
//                 //           << xpp[3] << "\n";
//                 // std::cout << "xmm: " << xmm[0] << "," << xmm[1] << "," << xmm[2] << ","
//                 //           << xmm[3] << "\n";

//                 const Float fact_nu_pp = (1.0 / 2.0) * eta(xpp, nu);
//                 const Float fact_nu_pm = (1.0 / 2.0) * eta(xpm, nu);
//                 const Float fact_nu_mp = (1.0 / 2.0) * eta(xmp, nu);
//                 const Float fact_nu_mm = (1.0 / 2.0) * eta(xmm, nu);

//                 // std::cout << "factors " << fact_nu_pp - fact_nu_mm << "\n";

//                 phi(x) += +fact_eta_x_mu * fact_nu_pm * (*U)(x, mu) *
//                           (*U)(xpm, nu).dagger() * psi(xpm);
//                 phi(x) +=
//                   -fact_eta_x_mu * fact_nu_pp * (*U)(x, mu) * (*U)(xp, nu) * psi(xpp);
//                 phi(x) += -fact_eta_x_mu * fact_nu_mm * (*U)(xm, mu).dagger() *
//                           (*U)(xmm, nu).dagger() * psi(xmm);
//                 phi(x) += +fact_eta_x_mu * fact_nu_mp * (*U)(xm, mu).dagger() *
//                           (*U)(xm, nu) * psi(xmp);

//                 xpp[nu]--; // x+mu again
//                 xpm[nu]++; // x+mu again
//                 xmp[nu]--; // x-mu again
//                 xmm[nu]++; // x-mu again
//               }
//               // adding the contributions from m*(G(x,y) + G^{\dagger}(x,y))
//               phi(x) += +m * fact_eta_x_mu * (*U)(x, mu) * psi(xp);
//               phi(x) += -m * fact_eta_x_mu * (*U)(xm, mu).dagger() * psi(xm);

//               phi(x) += +m * fact_eta_xm_mu * (*U)(xm, mu).dagger() * psi(xm);
//               phi(x) += -m * fact_eta_xp_mu * (*U)(x, mu) * psi(xp);

//               xm[mu]++; // =x again
//               xp[mu]--; // =x again

//               xpp[mu]--; // =x again
//               xpm[mu]--; // =x again
//               xmp[mu]++; // =x again
//               xmm[mu]++; // =x again
//             }
//             // adding the contribution from m^2*delta_{x,y}
//             phi(x) += std::pow(m, 2.0) * psi(x);
//           }
//         }
//       }
//     }

    return phi;
  }

  // returns the reslut of D*psi, where D is the Dirac operator
  template <class Float, class Type, class Group>
  spinor_lat<Float, Type>
  apply_D(gaugeconfig<Group> *U, const Float &m, const spinor_lat<Float, Type> &psi) {
    const size_t Lt = U->getLt(), Lx = U->getLx(), Ly = U->getLy(), Lz = U->getLz();
    const size_t nd = U->getndims();
    const std::vector<size_t> dims = psi.get_dims(); // vector of dimensions

    const int N = psi.size();
    spinor_lat<Float, Type> phi(dims, N);
#pragma omp parallel for
    for (int x0 = 0; x0 < Lt; x0++) {
      for (int x1 = 0; x1 < Lx; x1++) {
        for (int x2 = 0; x2 < Ly; x2++) {
          for (int x3 = 0; x3 < Lz; x3++) {
            const std::vector<int> x = {x0, x1, x2, x3};
            std::vector<int> xm = x, xp = x;
            for (size_t mu = 0; mu < nd; mu++) {
              const Float eta_x_mu = eta(x, mu);
              xm[mu]--; // x - mu
              xp[mu]++; // x + mu

              phi(x) += +(1.0 / 2.0) * eta_x_mu * (*U)(x, mu) * psi(xp);
              phi(x) += -(1.0 / 2.0) * eta_x_mu * (*U)(xm, mu).dagger() * psi(xm);

              xm[mu]++; // =x again
              xp[mu]--; // =x again
            }
            phi(x) += m * psi(x);
          }
        }
      }
    }
    return phi;
  }

  // returns the reslut of D^{\dagger}*psi, where D is the Dirac operator
  template <class Float, class Type, class Group>
  spinor_lat<Float, Type>
  apply_Ddag(gaugeconfig<Group> *U, const Float &m, const spinor_lat<Float, Type> &psi) {
    const size_t Lt = U->getLt(), Lx = U->getLx(), Ly = U->getLy(), Lz = U->getLz();
    const size_t nd = U->getndims();
    const std::vector<size_t> dims = psi.get_dims(); // vector of dimensions

    const int N = psi.size();
    spinor_lat<Float, Type> phi(dims, N);
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

              const Float eta_xm_mu = eta(xm, mu);
              const Float eta_xp_mu = eta(xp, mu);

              phi(x) += +(1.0 / 2.0) * eta_xm_mu * (*U)(xm, mu).dagger() * psi(xm);
              phi(x) += -(1.0 / 2.0) * eta_xp_mu * (*U)(x, mu) * psi(xp);

              xm[mu]++; // =x again
              xp[mu]--; // =x again
            }
            phi(x) += m * psi(x);
          }
        }
      }
    }
    return phi;
  }

} // namespace staggered
