// propagator.hpp
/**
 * @brief Dirac propagator routines
 * Computation of the Dirac propagator
 */

#pragma once
#ifndef Genz
#include "./staggered.hpp"

namespace staggered {

  /**
   * @brief pion correlator at time t and at rest (\vec{p}=\vec{0})
   * returns:
   * \sum_{\vec{x}} |D^{-1}(x|0)|^2, where x=(t, \vec{x})
   * (see eq. 6.32 of https://link.springer.com/book/10.1007/978-3-642-01850-3)
   * In order to use CG inversion algorithms (need for hermitian matrices),
   * we find the propagator S (inverse of the Dirac operator) as:
   * S = D^{-1} = D^{\dagger}*(D*D^{\dagger})^{-1}
   * @param U gauge configuration
   * @param m mass
   * @return std::vector vector of values for each t
   */
  template <class Group, class Float = double, class Type = std::complex<double>>
  std::vector<Float> C_pion(const gaugeconfig<Group> &U,
                            const Float &m,
                            const std::string &solver,
                            const Float &tol,
                            const size_t &verb,
                            const size_t &seed) {
    const size_t Lt = U.getLt(), Lx = U.getLx(), Ly = U.getLy(), Lz = U.getLz();
    const size_t nd = U.getndims();
    const nd_max_arr<size_t> dims = {Lt, Lx, Ly, Lz};

    std::vector<Float> C(Lt, 0.0); // correlator

    spinor_lat<Float, Type> source(dims, 0.0);
    source({0, 0, 0, 0}) = 1.0;
    const staggered::DDdag_matrix_lat<Float, Complex, Group> DDdag(U, m);

    const double Vs = (U.getVolume() / U.getLt()); // spatial volume
#pragma omp parallel for
    for (int x0 = 0; x0 < Lt; x0++) {
      if (x0 == 0) {
        continue; // contact divergence at x0==0 --> no need to compute it
      }
      for (int x1 = 0; x1 < Lx; x1++) {
        for (int x2 = 0; x2 < Ly; x2++) {
          for (int x3 = 0; x3 < Lz; x3++) {
            const nd_max_arr<int> x = {x0, x1, x2, x3};

            const spinor_lat<Float, Type> R = DDdag.inv(source, solver, tol, verb, seed);

            const Type res = apply_Ddag(U, m, R)(x);
            C[x0] += retrace(conj(res) * res);
          }
        }
      }
      C[x0] /= Vs;
    }
    return C;
  }

} // namespace staggered
#endif