// BiCGStab.hpp
/*
Simone Romiti - simone.romiti.1994@gmail.com

Bi Conjugate Gradient method
Reference: sec. 8.8.4 of
Degrand-Lattice Methods for Quantum Chromodynamics
https://www.worldscientific.com/worldscibooks/10.1142/6065
*/

#pragma once

#include <iostream>
#include <vector>

#include "CG_solver.hpp"

namespace BiCGStab {
  using namespace CG_solver;

  /**
   * Conjugate Gradient method class
   * @Float = type for the residual
   * @T = type stored inside vectors and matrices.
   * The method T::real(), which gives the real part, must have been implemented
   * @LAmatrix, @LAvector = types of matrix and vector
   */
  template <class Float, class T, class LAmatrix, class LAvector>
  class LinearBiCGStab : public LinearCG_solver<Float, T, LAmatrix, LAvector> {
  public:
    LinearBiCGStab() {}
    ~LinearBiCGStab() {}

    LinearBiCGStab(const LAmatrix &_A, const LAvector &_b) { this->init_system(_A, _b); }

    void
    solve(const LAvector &x0, const Float &tol = 1e-15, const size_t &verbosity = 0) {
      this->check_initial_guess(x0);

      const LAvector b = (*this).b;
      const LAmatrix A = (*this).A;
      LAvector xk = x0; // initial vector
      LAvector r0 = b - A * xk; // residual vector
      LAvector rk = r0;
      LAvector pk = rk;
      Float rk_norm = rk.norm();

      int num_iter = 0;
      ((*this).curve_x).push_back(xk);
      while (rk_norm > tol) {
        const LAvector apk = A * pk;
        const T r0rk = r0.dot(rk);

        const T alpha = r0.dot(rk) / r0.dot(apk);
        const LAvector s = rk - alpha * apk;
        const LAvector t = A * s;

        const T omega = t.dot(s) / t.dot(t);

        rk = s - omega * (A * s);
        xk = xk + omega * s + alpha * pk;

        const T beta = (r0.dot(rk) / r0rk) * (alpha / omega);
        pk = rk + beta * pk - beta * omega * apk;

        num_iter += 1;
        ((*this).curve_x).push_back(xk);
        rk_norm = rk.norm();

        if (verbosity > 1) {
          std::cout << "Iteration: " << num_iter << "\n";
          if (verbosity > 2) {
            std::cout << "x = ";
            print_LAvector<T>(xk, ",");
          }
          std::cout << "residual rk.norm() = " << rk_norm << "\t;\t";
          std::cout << "(b - A*x).norm() = " << (b - A * xk).norm() << "\n";
        }
      }

      const Float ex_res = (b - A * xk).norm(); // exact residual

      if (rk_norm < ex_res) {
        if (verbosity > 1) {
          std::cout << "\nRoundoff error detected: (rk - A*pk) != (b - A*xk).norm() . "
                       "Repeating the  the BiCGStab.\n\n";
        }

        // LinearBiCGStab<Float, T, LAmatrix, LAvector> Lloop(A, b-A*xk);
        // LAvector dx0;
        // dx0.resize(xk.size());
        // Lloop.solve(dx0, tol, verbosity);
        // xk = xk + Lloop.get_solution();
        // std::cout << "ciao2\n";
        this->solve(xk, tol, verbosity); // repeat the CG starting from the found solution
                                         // until (rk - apk) == ((A*xk) - b).norm()
      } else {
        if (verbosity > 0) { // printing the solution
          std::cout << "\nSolution: \t x = ";
          print_LAvector<T>(xk, ",");
        }

        (*this).x = xk; // storing the solution
        (*this).solved = true; // system has been solved
      }
    }

  }; // class LinearBiCGStab

} // namespace BiCGStab
