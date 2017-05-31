#include "gaugeconfig.hh"
#include "md_update.hh"
#include "compute_hamiltonian.hh"
#include<vector>
#include<random>

using std::vector;

template<class URNG> bool md_update(gaugeconfig &U,
                                    URNG const &engine, 
                                    md_params const &params) {
  size_t Volume = U.getVolume();
  vector<double> momenta;
  momenta.resize(Volume*4*3);
  std::normal_distribution<double> normal(0., 1.);
  std::uniform_real_distribution<double> uniform(0., 1.);
  for(int i = 0; i < Volume*4*3; i++) {
    momenta[i] = normal(engine);
  }

  const double H_old = compute_hamiltonian(momenta, U);
  gaugeconfig U_old(U);
  bool accepted = true;
  
  double delta_H = H_old - compute_hamiltonian(momenta, U);
  if(delta_H > 0) {
    if(uniform(engine) > exp(-delta_H)) {
      accepted = false;
    }
  }
  if(!accepted) {
    U = U_old;
  }

  return accepted;
}
