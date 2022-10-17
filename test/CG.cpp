// CG.cpp

#include <complex>

#include "CG.hpp"
#include "LA.hpp"
// #include "./cg_invert.hpp" // DEPRECATED

typedef std::complex<double> Type;

const Type i(0.0, 1.0);

typedef LA::LAmatrix<Type> LAmatrix;
typedef LA::LAvector<double, Type> LAvector;

int main(int argc, char const *argv[]) {

  std::cout << "running CG.cpp\n"; 

  LAmatrix A(0, 0);
  A.add_row((std::vector<Type>){1.0, 2.0+i, 1.0});
  A.add_row((std::vector<Type>){2.0-i, 5.0, -6.0 + 6.0*i});
  A.add_row((std::vector<Type>){1.0, -6.0 - 6.0*i, -6.0});

//   A.add_row((std::vector<Type>){1.0, 2.0, 1.0});
//   A.add_row((std::vector<Type>){2.0, 5.0, -6.0});
//   A.add_row((std::vector<Type>){1.0, -6.0, -6.0});

  std::cout << "A=" << '\n';
  CG::print_LAmatrix<Type, LAmatrix, LAvector>(A, ",");
  const LAvector x_star = (std::vector<Type>){10.0, 11.0, -1.0};
  const LAvector x0 = (std::vector<Type>){5e+5, -1e+7, 0.0};

  const LAvector b = A * x_star;

  std::cout << "The solution found by the cg is:" << '\n';

  CG::LinearCG<double, Type, LAmatrix, LAvector> LCG(A, b);
  LCG.solve(x0, 1e-16, 1);

  std::cout << "Exact solution: x = ";
  CG::print_LAvector<Type, LAvector>(x_star, ",");

  // // std::cin.get();
  // std::cout << "Computing the inverse:" << '\n';
  // const int n_hits = 15;
  // CG::inverter<Type, LAmatrix, LAvector> CGI(A, n_hits, 1);
  // CGI.invert(1e-15, 1);
  // std::cout << "The solution found with the inverter is:" << '\n';
  // LAmatrix B = CGI.get_inverse();
  // CG::print_LAvector<Type, LAvector>(B*b, ",");

  return 0;
}
