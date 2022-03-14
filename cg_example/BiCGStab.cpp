// cg.cpp

#include <complex>

#include "../BiCGStab.hpp"
#include "./LA.hpp"
// #include "./cg_invert.hpp" // DEPRECATED

typedef std::complex<double> Type;

const Type i(0.0, 1.0);

typedef LA::LAmatrix<Type> LAmatrix;
typedef LA::LAvector<double, Type> LAvector;

int main(int argc, char const *argv[]) {
  LAmatrix A(0, 0);
  A.add_row((std::vector<Type>){1.0, 2.0+i, 1.0});
  A.add_row((std::vector<Type>){2.0-i, 5.0, -6.0 + 6.0*i});
  A.add_row((std::vector<Type>){1.0, -6.0 - 6.0*i, -6.0});

//   A.add_row((std::vector<Type>){1.0, 2.0, 1.0});
//   A.add_row((std::vector<Type>){2.0, 5.0, -6.0});
//   A.add_row((std::vector<Type>){1.0, -6.0, -6.0});

  std::cout << "A=" << '\n';
  BiCGStab::print_LAmatrix<Type, LAmatrix, LAvector>(A, ",");
  const LAvector x_star = (std::vector<Type>){10.0, 11.0, -1.0};
  const LAvector x0 = (std::vector<Type>){5e+5, -1e+7, 0.0};

  const LAvector b = A * x_star;

  std::cout << "The solution found by the BiCGStab is:" << '\n';

  BiCGStab::LinearBiCGStab<double, Type, LAmatrix, LAvector> LBiCGStab(A, b);
  LBiCGStab.solve(x0, 1e-25, 2);

  std::cout << "Exact solution: x = ";
  BiCGStab::print_LAvector<Type, LAvector>(x_star, ",");


  return 0;
}
