// links.cpp

// #include "fermions/cg.hpp"
#include "fermions/staggered.hpp"
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char const *argv[]) {
  const int Lx = 4, Ly = 4, Lz = 4, Lt = 4;
  const int ndims = 4;
  const double beta = 1.5;
  gaugeconfig<u1> U(Lx, Ly, Lz, Lt, ndims, beta);
  hotstart<u1>(U, 187238, true);

  std::pair<std::complex<double> *, size_t *> p =
    staggered::dirac_op_ptr<std::complex<double>>(U, 0.1);
  std::complex<double> *D = p.first;
  size_t *col_idx = p.second;

  // Create a sparse matrix
  const size_t N_c = 1;
  const size_t N = U.getSize() * N_c;
  Eigen::SparseMatrix<std::complex<double>> M(N, N);
  for (size_t i = 0; i < N; i++) {
    M.insert(3 * i, col_idx[3 * i]) = D[3 * i];
    M.insert(3 * i, col_idx[3 * i + 1]) = D[3 * i + 1];
    M.insert(3 * i, col_idx[3 * i + 2]) = D[3 * i + 2];
  }

  // Make sure to finalize the matrix after inserting coefficients
  M.finalize();
  M = M*(M.adjoint());

//   for (size_t i = 0; i < N; i++) {
//     for (size_t j = 0; j < N; j++) {
//       std::cout << M(i, j) << " ";
//     }
//     std::cout << "\n";
//   }

  Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>, Eigen::COLAMDOrdering<int>>
    solver;
  // fill A and b;
  // Compute the ordering permutation vector from the structural pattern of A
  solver.analyzePattern(M);
  // Compute the numerical factorization
  solver.factorize(M);

  // Compute the determinant
  std::cout << "determinant: \n";
  std::cout << solver.absDeterminant() << "\n";

  return 0;
}
