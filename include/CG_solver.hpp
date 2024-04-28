// CG_solver.hpp
/*
Simone Romiti - simone.romiti.1994@gmail.com

This file defines a general purpose, templated class that is the basis for conjugate
gradient (CG) algorithms for solving linear systems. Please recall that convergence is
guaranteed only for invertible hermitean matrices.

It is assumed there are types LAvector and
LAmatrix (LA=Linear Algebra) such that are defined:
  1.  template<class T> LAvector operator *(const LAmatrix&, const LAvector&) // A*v
  2.  template<class T> LAvector operator -(const LAvector&, const LAvector&) // a-b
  3.  template<class T> LAvector operator +(const LAvector&, const LAvector&) // a+b
  4.  template<class T> LAvector operator -(const LAvector&) // -v
  5.  template<class T> LAvector operator *(const T& lambda, const LAvector& v) //
lambda*v
  6.  template<class T> T operator *(const LAvector& a, const LAvector& b) //
a_i^{\dagger} * b_i : complex dot product
  7.  template<class Float, class LAvector> norm(const LAvector& v)
  8.  template<class Float, class T, class LAmatrix, LAvector> size_t LAmatrix::rows()
const; // numer of rows
  9.  template<class Float, class T, class LAmatrix, LAvector> size_t LAmatrix::cols()
const; // numer of columns
  10. template<class Float, class T, class LAmatrix, LAvector> size_t LAvector::size()
const; // vector's length
  11. template<class Float, class T, class LAmatrix, LAvector> size_t LAvector::dot(const
LAvector& w) const; // complex dot product <v, w>
  12. template<class Float, class T, class LAmatrix, LAvector> void LAmatrix::operator
=(const LAmatrix&); // add row to the matrix This kind of implementation allows for
user-defined parallelization of matrix-vector multiplication, norm evaluation, etc.
*/

#pragma once

#include <iostream>
#include <vector>

namespace CG_solver {

  template <class T, class LAvector>
  void print_LAvector(const LAvector &v, const std::string &sep = ",") {
    const int N = v.size();
    std::cout << "{";
    for (int i = 0; i < N; ++i) {
      std::cout << v[i] << sep << " ";
    }
    std::cout << "}\n";
  }

  template <class T, class LAmatrix, class LAvector>
  void print_LAmatrix(const LAmatrix &m, const std::string &sep = ",") {
    std::cout << "{";
    const int nr = m.rows();
    for (int i = 0; i < nr; ++i) {
      print_LAvector<T, LAvector>(m[i], sep);
    }
    std::cout << "}\n";
  }

  /**
   * Base class for Conjugate Gradient methods
   * @Float = type for the residual
   * @T = type stored inside vectors and matrices.
   * The method T::real(), which gives the real part, must have been implemented
   * @LAmatrix, @LAvector = types of matrix and vector
   */
  template <class Float, class T, class LAmatrix, class LAvector> class LinearCG_solver {
  protected:
    LAmatrix A;
    LAvector b;
    LAvector x; // solution to A*x = b
    bool sys_init = false; // true when A and b are initialized
    std::vector<LAvector> curve_x; // trajectory
    bool solved = false; // true when x has been found

    void check_solved() const {
      if (!solved) {
        std::cerr << "Error. Can't get the trajectory: solution hasn't been found yet. "
                     "Aborting.\n";
       std::abort();
      }
    }

  public:
    LinearCG_solver() {}
    ~LinearCG_solver() {}

    void init_system(const LAmatrix &_A, const LAvector &_b) {
      if (_A.rows() != _b.size()) {
        std::cout
          << "Error. Invalid linear system. Check matrix and vector sizes. Aborting. \n";
       std::abort();
      }

      A = _A;
      b = _b;
      sys_init = true;
    }

    LinearCG_solver(const LAmatrix &_A, const LAvector &_b) { this->sys_init(_A, _b); }

    void check_initial_guess(const LAvector &x0) const {
      if (A.cols() != x0.size()) {
        std::cerr
          << "Error. Invalid starting vector x0. Check the number of components:\n";
        std::cerr << "A.cols(): " << A.cols() << "\n";
        std::cerr << "x0.cols(): " << x0.size() << "\n";
        std::cerr << "Aborting.\n";
       std::abort();
      }
    }

    LAmatrix get_trajectory() {
      this->check_solved();
      return curve_x;
    }

    LAvector get_solution() {
      this->check_solved();
      return x;
    }

  }; // class LinearCG_solver

} // namespace CG_solver
