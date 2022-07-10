// links.cpp

#include <complex>

#include "../include/adjointfield.hh"
#include "../include/loops.hpp"
#include "../include/operators.hpp"
#include "../include/u1.hh"

template <class T> using nd_max_arr = spacetime_lattice::nd_max_arr<T>;

int main(int argc, char const *argv[]) {
  size_t d = 2, l = 4;

  std::cout << d << " " << l << "\n";

  links::closed_paths P(d, 0, l);
  std::cout << P.size() << "\n";

  // std::cout << "Printing loop:\n";
  // P.print();
  // std::cout << "Aborting\n";
  // std::abort();

  const int nd = 2;
  gaugeconfig<_u1> U(4, 1, 1, 8, nd, 1.2345678);
  hotstart(U, 1391693, true); // random gauge config

  //  nd_max_arr<int> x = {0,0,0,0}; // origin
  double P1 = 0.0;
  nd_max_arr<int> xmm = {-1, -1, 0, 0}, xm0 = {-1, 0, 0, 0}, xmp = {-1, 1, 0, 0};
  nd_max_arr<int> x0m = {0, -1, 0, 0}, x00 = {0, 0, 0, 0}, x0p = {0, 1, 0, 0};
  nd_max_arr<int> xpm = {+1, -1, 0, 0}, xp0 = {+1, 0, 0, 0}, xpp = {+1, 1, 0, 0};
  P1 += retrace(U(x00, 0) * U(xp0, 1) * U(x0p, 0).dagger() * U(x00, 1).dagger());
  P1 += retrace(U(xm0, 0) * U(x00, 1) * U(xmp, 0).dagger() * U(xm0, 1).dagger());
  P1 += retrace(U(xmm, 0) * U(x0m, 1) * U(xm0, 0).dagger() * U(xmm, 1).dagger());
  P1 += retrace(U(x0m, 0) * U(xpm, 1) * U(x00, 0).dagger() * U(x0m, 1).dagger());
  std::cout << "P1 : " << P1 << "\n";

  // const double dims_fact = spacetime_lattice::num_pLloops_half(U.getndims() - 1);
  links::links_loops<_u1> LL(U, x00, 0, 4);
  LL.find_all_shorter();
  std::vector<_u1> loops = LL.get_loops();

  double P2 = 0.0;
  for (size_t i = 0; i < loops.size(); i++) {
    P2 += retrace(loops[i]) / 2;
  }
  std::cout << "P2 : " << P2 << "\n";

  if (P1 == P2) {
    std::cout << "Success: P1==P2\n";
  }

  return 0;
}
