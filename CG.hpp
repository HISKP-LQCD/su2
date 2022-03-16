// CG.hpp
/*
Simone Romiti - simone.romiti@uni-bonn.de

Conjugate Gradient method
Reference: sec. 8.8.2 of
Degrand-Lattice Methods for Quantum Chromodynamics
https://www.worldscientific.com/worldscibooks/10.1142/6065
*/

#pragma once

#include <iostream>
#include <vector>

#include "CG_solver.hpp"

namespace CG {
  using namespace CG_solver;

  /**
   * Conjugate Gradient method class
   * @Float = type for the residual
   * @T = type stored inside vectors and matrices.
   * The method T::real(), which gives the real part, must have been implemented
   * @LAmatrix, @LAvector = types of matrix and vector
   */
  template <class Float, class T, class LAmatrix, class LAvector>
  class LinearCG : public LinearCG_solver<Float, T, LAmatrix, LAvector> {
  public:
    LinearCG() {}
    ~LinearCG() {}

    LinearCG(const LAmatrix &_A, const LAvector &_b) { this->init_system(_A, _b); }

    void
    solve(const LAvector &x0, const Float &tol = 1e-15, const size_t &verbosity = 0) {
      this->check_initial_guess(x0);

      const LAvector b = (*this).b;
      const LAmatrix A = (*this).A;

      LAvector xk = x0; // initial vector
      LAvector rk = b - A * xk; //(A * xk) - b; // residual vector
      LAvector pk = rk;
      Float rk_norm = rk.norm();

      int num_iter = 0;
      (*this).curve_x.push_back(xk);
      while (rk_norm > tol) {
        const LAvector apk = A * pk;
        const T rkrk = rk.dot(rk);

        const T alpha = rk.dot(rk) / pk.dot(apk);
        xk = xk + (alpha * pk);
        rk = rk - (alpha * apk);

        const T beta = rk.dot(rk) / rkrk;
        pk = rk + (beta * pk);

        num_iter += 1;
        (*this).curve_x.push_back(xk);
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
                       "Repeating the  the CG.\n\n";
        }
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

  }; // class LinearCG

} // namespace CG
