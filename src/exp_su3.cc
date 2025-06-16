/**
 * @file exp_su3.cc
 * @author Simone Romiti (simone.romiti.1994@gmail.com)
 * @brief exponentiation of an element of the su(3) algebra: A=a_k*sigma_k -->
 * exp(i*a_k*lambda_k)
 * @version 0.1
 * @date 2023-02-21
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "adjointfield.hh"
#include "exp_gauge.hh"
#include "su3.hh"

#include <cmath>
#include <complex>
#include <vector>

/**
 * @brief exp(i*a_k*lambda_k), where the lambda_k are the Gell-Mann matrices and a_k are
 * real numbers.
 *
 * eq. (4) of https://arxiv.org/pdf/1508.00868v2.pdf
 * the Gell-Mann decomposition coefficients have been computed with sympy
 *
 * @param x
 * @return _su2
 */
_su3 exp(const adjointsu3<double> &x) {
  const double sqrt_of_3 = sqrt(3.0); // avoid computing it many times

  std::array<double, 8> arr = x.get_arr();
  double theta = 0.0;
  for (size_t i = 0; i < 8; i++) {
    theta += (arr[i] * arr[i]);
  }
  theta = std::sqrt(theta);

  // normalised vector
  std::array<double, 8> n = arr;
  for (size_t i = 0; i < 8; i++) {
    n[i] /= theta;
  }

  // determinant of H
  const double detH =
    2.0 * sqrt_of_3 * (n[0] * n[0]) * n[7] / 3.0 + 2.0 * n[0] * n[3] * n[5] +
    2.0 * n[0] * n[4] * n[6] + 2.0 * sqrt_of_3 * (n[1] * n[1]) * n[7] / 3.0 -
    2.0 * n[1] * n[3] * n[6] + 2.0 * n[1] * n[4] * n[5] +
    2.0 * sqrt_of_3 * (n[2] * n[2]) * n[7] / 3.0 + n[2] * (n[3] * n[3]) +
    n[2] * (n[4] * n[4]) - n[2] * (n[5] * n[5]) - n[2] * (n[6] * n[6]) -
    sqrt_of_3 * (n[3] * n[3]) * n[7] / 3.0 - sqrt_of_3 * (n[4] * n[4]) * n[7] / 3.0 -
    sqrt_of_3 * (n[5] * n[5]) * n[7] / 3.0 - sqrt_of_3 * (n[6] * n[6]) * n[7] / 3.0 -
    2.0 * sqrt_of_3 * (n[7] * n[7] * n[7]) / 9.0;
  const double phi = (1.0 / 3.0) * (acos((3.0 / 2.0) * sqrt_of_3 * detH) - M_PI / 2.0);

  const Complex I(0.0, 1.0); // imaginary unit

  // 1st and 2nd row of H
  const std::array<Complex, 3> H_1 = {n[2] + sqrt_of_3 * n[7] / 3.0, n[0] - I * n[1],
                                      n[3] - I * n[4]};
  const std::array<Complex, 3> H_2 = {n[0] + I * n[1], -n[2] + sqrt_of_3 * n[7] / 3.0,
                                      n[5] - I * n[6]};

  // 1st and 2nd row of H^2
  const std::array<Complex, 3> H2_1 = {
    (n[0] - I * n[1]) * (n[0] + I * n[1]) + std::pow(n[2] + sqrt_of_3 * n[7] / 3.0, 2.0) +
      (n[3] - I * n[4]) * (n[3] + I * n[4]),
    (n[0] - I * n[1]) * (-n[2] + sqrt_of_3 * n[7] / 3.0) +
      (n[0] - I * n[1]) * (n[2] + sqrt_of_3 * n[7] / 3.0) +
      (n[3] - I * n[4]) * (n[5] + I * n[6]),
    -2 * sqrt_of_3 * n[7] * (n[3] - I * n[4]) / 3.0 +
      (n[0] - I * n[1]) * (n[5] - I * n[6]) +
      (n[2] + sqrt_of_3 * n[7] / 3.0) * (n[3] - I * n[4])};
  const std::array<Complex, 3> H2_2 = {
    (n[0] + I * n[1]) * (-n[2] + sqrt_of_3 * n[7] / 3.0) +
      (n[0] + I * n[1]) * (n[2] + sqrt_of_3 * n[7] / 3.0) +
      (n[3] + I * n[4]) * (n[5] - I * n[6]),
    (n[0] - I * n[1]) * (n[0] + I * n[1]) +
      std::pow(-n[2] + sqrt_of_3 * n[7] / 3.0, 2.0) +
      (n[5] - I * n[6]) * (n[5] + I * n[6]),
    -2 * sqrt_of_3 * n[7] * (n[5] - I * n[6]) / 3.0 +
      (n[0] + I * n[1]) * (n[3] - I * n[4]) +
      (-n[2] + sqrt_of_3 * n[7] / 3.0) * (n[5] - I * n[6])};

  std::array<Complex, 3> u = {0.0, 0.0, 0.0}, v = {0.0, 0.0, 0.0};
  for (size_t k = 0; k <= 2; k++) {
    const double pk = phi + (2.0 / 3.0) * M_PI * k;
    const double arg_num_k = (2.0 / sqrt_of_3) * theta * sin(pk);
    const Complex num_k = std::exp(I * arg_num_k);
    const double den_k = 1 - 2.0 * cos(2.0 * pk);
    const Complex fact_k = num_k / den_k;
    for (size_t i = 0; i < 3; i++) {
      const Complex u_ik = H2_1[i] + (2.0 / sqrt_of_3) * H_1[i] * sin(pk);
      u[i] += u_ik * fact_k;

      const Complex v_ik = H2_2[i] + (2.0 / sqrt_of_3) * H_2[i] * sin(pk);

      v[i] += v_ik * fact_k;
    }
    // contibution from the identity
    u[0] += -(1.0 / 3.0) * (1.0 + 2.0 * cos(2.0 * pk)) * fact_k;
    v[1] += -(1.0 / 3.0) * (1.0 + 2.0 * cos(2.0 * pk)) * fact_k;
  }

  _su3 res(u, v);
  return res;
}
