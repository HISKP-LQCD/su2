// dirac.cpp

#include <complex>
#include <fstream>
#include <iostream>
#include <string>

// #include <Eigen/LU>
// #include <Eigen/Sparse>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "fermions/staggered.hpp"

namespace ublas = boost::numeric::ublas;

int main(int argc, char const *argv[]) {
  const int Lx = 2, Ly = 2, Lz = 1, Lt = 2;
  const int ndims = 3;
  const double beta = 1.5;
  std::cout << "Initializing the configuration\n";
  gaugeconfig<u1> U(Lx, Ly, Lz, Lt, ndims, beta);
  hotstart<u1>(U, 187238, true);

  std::cout << "Building the Dirac operator\n";
  std::pair<std::complex<double> *, size_t *> p =
    staggered::dirac_op_ptr<std::complex<double>>(U, 0.1);
  std::complex<double> *D = p.first;
  size_t *col_idx = p.second;

  // Create a sparse matrix
  const size_t N_c = 1;
  const size_t N = U.getSize() * N_c;

  std::cout << "Converting the Dirac operator to the ublas format\n";
  // ublas doesn't support LU factorization for sparse complex matrices
  // therefore I store each complex number a+ib as the matrix ((a,-b),(b,a))
  ublas::compressed_matrix<double> M(2 * N, 2 * N);
  for (size_t i = 0; i < N; i++) {
    for (size_t k = 0; k < 3; k++) {
      const size_t j = col_idx[3 * i + k];
      const std::complex<double> Mij = D[3 * i + k];
      const double reMij = std::real(Mij), imMij = std::imag(Mij);
      M(2 * i + 0, 2 * j + 0) = reMij;
      M(2 * i + 1, 2 * j + 1) = reMij;
      M(2 * i + 1, 2 * j + 0) = -imMij;
      M(2 * i + 0, 2 * j + 1) = imMij;
    }
  }

  std::cout << "LU factorization\n";
  ublas::compressed_matrix<double> LU = M;
  ublas::permutation_matrix<std::size_t> pm(M.size1());
  ublas::lu_factorize(LU, pm);

  std::cout << "Calculating determinant from LU decomposition\n";
  double determinant = 1.0;
  for (std::size_t i = 0; i < LU.size1(); ++i) {
    determinant *= (i == pm(i)) ? 1 : -1;
    determinant *= LU(i, i);
  }
  std::cout << "Determinant: " << determinant << std::endl;

  return 0;
}
