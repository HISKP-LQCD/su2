#include"gradient_flow.hh"
#include"su2.hh"
#include"gaugeconfig.hh"
#include"adjointfield.hh"
#include"gauge_energy.hh"
#include"energy_density.hh"
#include"monomial.hh"
#include"gaugemonomial.hh"
#include"hamiltonian_field.hh"
#include"update_gauge.hh"

#include<string>
#include<fstream>

void runge_kutta(hamiltonian_field<double> &h, monomial<double> &SW, const double eps) {

  double zfac[5] = { (-17.0)/(36.0), (8.0)/(9.0), (-3.0)/(4.0)};
  double expfac[3] = {-36.0/4./17.0, 1., -1.};

  // before the three steps zero the derivative field
  zeroadjointfield(*(h.momenta));

  for(int f = 0; f < 3; f++) {
    // add to *(h.momenta) 
    SW.derivative(*(h.momenta), h, zfac[f]/h.U->getBeta());
    // where does the '-' come from?
    // update the flowed gauge field Vt
    update_gauge(h, -eps*expfac[f]);
  }
  return;
}

void gradient_flow(gaugeconfig &U, std::string const &path) {
  double t[3], P[3], E[3];
  double eps = 0.01;
  std::ofstream os(path, std::ios::out);

  for(unsigned int i = 0; i < 3; i++) {
    t[i] = 0.;
    P[i] = 0.;
    E[i] = 0.;
  }
  P[2] = gauge_energy(U)/U.getVolume()/N_c/6.;
  E[2] = energy_density(U);

  gaugeconfig Vt(U);
  adjointfield<double> deriv(U.getLs(), U.getLt());
  hamiltonian_field<double> h(deriv, Vt);

  gaugemonomial<double> SW(0);

  while(t[1] < 9.99) {
    t[0] = t[2];
    P[0] = P[2];
    E[0] = E[2];
    for(unsigned int x0 = 1; x0 < 3; x0++) {
      t[x0] = t[x0-1] + eps;
      runge_kutta(h, SW, eps);
      P[x0] = gauge_energy(Vt)/U.getVolume()/N_c/6.;
      E[x0] = energy_density(Vt);
    }
    double tsqP = t[1]*t[1]*2*N_c*6.*(1-P[1]);
    double tsqE = t[1]*t[1]*E[1];
    os << t[1] << " "  << P[1] << " " << 2*N_c*6.*(1.-P[1]) << " " << tsqP << " " << E[1] << " " << tsqE << std::endl;
  }

  return;
}
