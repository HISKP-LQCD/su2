#include"compute_hamiltonian.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"su2.hh"
#include<vector>
#include<complex>

double compute_hamiltonian(std::vector<double> const &momenta,
                           gaugeconfig &U) {
  double H = 0.;
  for(size_t i = 0; i < U.getVolume()*4*3; i++) {
    H += 0.5*momenta[i]*momenta[i];
  }
  H += U.getBeta()/N_c*gauge_energy(U);
  return(H);
}
