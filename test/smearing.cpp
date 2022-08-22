// cg.cpp

#include <complex>
#include <iostream>

#include "geometry.hh"
#include "random_gauge_trafo.hh"

#include "smearape.hh"
template <class T> using nd_max_arr = spacetime_lattice::nd_max_arr<T>;

// typedef std::complex<double> Type;

bool g_heat = true; // hot configuration to start from
const size_t L = 8, T = 16; // lattice sizes
const size_t ndims = 3;

const nd_max_arr<int> x = {5, 2, 0, 0};
const size_t mu = 1;

void print_phase(const _u1 &u, const std::string &prep = "") {
  const std::complex<double> t1 = trace(u);
  if (std::abs(t1) != 1.0) {
    std::cerr << "Error. Something's wrong. You didn't pass a pure phase. absolute value "
                 "!= 1.0.\n";
    std::cerr << "Aborting.";
    std::abort();
  }
  std::cout << prep << t1 << " " << std::arg(t1) << "\n";
}

int main(int argc, char const *argv[]) {
  std::cout.precision(15);
  std::cout << "running smearing.cpp\n\n";

  // random gauge configuration
  gaugeconfig<_u1> U0(L, L, L, T, ndims, 2.1);
  hotstart(U0, 52125, g_heat);

  const gaugeconfig<_u1> U = U0;

  std::cout << "\t\t printing format:  \n";
  std::cout << "\t\t-------------------\n";
  std::cout << "\t\te^{i*alpha} | alpha\n";
  std::cout << "\t\t-------------------\n";
  std::cout << "\n";

  print_phase(U(x, mu), "Initial configuration   : ");

  gaugeconfig<_u1> U1 = U; // copy of the initial configuration

  const double alpha1 = 1.0; // no smearing
  spatial_APEsmearing_u1<double>(U1, alpha1);
  print_phase(U1(x, mu), "no smearing (alpha=1.0) : ");

  _u1 Udag_U1 = U(x, mu).dagger() * U1(x, mu);
  print_phase(Udag_U1, "Udag*U1 : ");

  gaugeconfig<_u1> U2 = U; // copy of the initial configuration
  const double alpha2 = 0.0; // full smearing
  spatial_APEsmearing_u1<double>(U2, alpha2);
  print_phase(U2(x, mu), "full smearing (alpha=0.0) : ");
  std::cout << "\n";

  std::cout << "Udag*U2\n";
  _u1 Udag_U2 = U(x, mu).dagger() * U2(x, mu);
  print_phase(Udag_U2);
  std::cout << "\n";

  std::cout << "Applying a random gauge transformation many times:\n";
  std::cout << "Udag*U2 should not change. It is a sum of plaquettes.\n";

  gaugeconfig<_u1> U_rnd = U; // copy of the initial configuration
  gaugeconfig<_u1> U2_rnd = U2; // copy of the initial smeared configuration

  for (size_t i = 0; i < 3; i++) {
    random_gauge_trafo(U_rnd, 654321);
    random_gauge_trafo(U2_rnd, 654321);
    Udag_U2 = U_rnd(x, mu).dagger() * U2_rnd(x, mu);
    std::cout << "i " << i << "\n";
    print_phase(U2_rnd(x, mu), "U2(x, mu) : ");
    print_phase(Udag_U2, "Udag*U2 : ");
  }

  return 0;
}
