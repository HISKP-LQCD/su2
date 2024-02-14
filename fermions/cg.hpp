// routines for the Conjugate-Gradient method

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/LU> 
#include <Eigen/Sparse>

namespace fermions {

  template <class Float>
  void cg_solve_eigen(const Eigen::SparseMatrix<Float> &A,
                      const Eigen::VectorXd &b,
                      const int &MaxIterations,
                      const double &tolerance) {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<Float>> cg;
    cg.setMaxIterations(MaxIterations);
    cg.setTolerance(tolerance);
    cg.compute(A);

    Eigen::VectorXd x(A.rows());
    x.setZero();
    x = cg.solve(b);
    if (cg.info() != Eigen::Success) {
      std::cerr << "Conjugate Gradient solver failed to converge!" << std::endl;
      std::abort();
    }
    std::cout << x << std::endl;
  }
} // namespace fermions