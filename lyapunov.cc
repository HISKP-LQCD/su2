//

template<class T> void compute_lyapunov(gaugeconfig &U, std::string const &path, integrator<T> &md_integ, const int d = 12) {
  adjointfield<T> momenta(U.getLs(), U.getLt());
  // generate standard normal distributed random momenta
  // normal distribution checked!
  momenta = initnormal<URNG, T>(engine, U.getLs(), U.getLt());
  adjointfield<T> momenta2(momenta);

  gaugeconfig U2(U.getLs(), U.getLt(), U.getbeta());

  return;
}
