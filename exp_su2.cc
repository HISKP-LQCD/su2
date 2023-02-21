/**
 * @file exp_su2.cc
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief exponentiation of an element of the su(2) algebra: A=alpha_k*sigma_k -->
 * exp(i*alpha_k*sigma_k)
 * @version 0.1
 * @date 2023-02-21
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <cmath>
#include <complex>
#include <vector>

#include "adjointfield.hh"
#include "exp_gauge.hh"
#include "su2.hh"

/**
 * @brief exp(i*alpha_k*sigma_k), where the sigma_k are the Pauli matrices and alpha_k are
 * real numbers.
 *
 * The code uses the explicit formula following by (4.23) of Gattringer&Lang
 * https://link.springer.com/book/10.1007/978-3-642-01850-3
 *
 * @param x
 * @return _su2
 */
_su2 exp(adjointsu2<double> const &x) {
  double a = x.geta(), b = x.getb(), c = x.getc();
  // normalise the vector
  const double alpha = sqrt(a * a + b * b + c * c);
  std::vector<double> n = {a / alpha, b / alpha, c / alpha};

  const double salpha = sin(alpha);
  _su2 res(Complex(cos(alpha), salpha * n[2]), salpha * Complex(n[1], n[0]));
  return res;
}
