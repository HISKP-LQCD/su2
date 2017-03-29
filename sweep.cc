#include"gaugeconfig.hh"
#include"random_su2.hh"
#include"get_staples.hh"
#include"sweep.hh"
#include<random>
#include<vector>

double sweep(gaugeconfig &U, const size_t seed, const double delta, const size_t N_hit, const double beta) {
  std::mt19937 engine(seed);
  std::uniform_real_distribution<double> uniform(0., 1.);

  size_t rate = 0;
  su2 rU(0., 0.);
  std::vector<size_t> x = {0, 0, 0, 0};
  for(x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for(x[1] = 0; x[1] < U.getLs(); x[1]++) {
      for(x[2] = 0; x[2] < U.getLs(); x[2]++) {
        for(x[3] = 0; x[3] < U.getLs(); x[3]++) {
          for(size_t mu = 0; mu < 4; mu++) {
            su2 K = get_staples(U, x, mu);
            for(size_t n = 0; n < N_hit; n++) {
              bool accept = false;
              random_su2(rU, engine, delta);
              double deltaS = beta/static_cast<double>(N_c)*(trace<su2>(U(x, mu) * K) - trace(U(x, mu) * rU * K));
              if(deltaS < 0) accept = true;
              else accept = (uniform(engine) < exp(-deltaS));
              if(accept) {
                U(x, mu) = U(x, mu) * rU;
                U(x, mu).rescale();
                rate += 1;
              }
            }
          }
        }
      }
    }
  }
  return( double(rate)/double(N_hit)/double(U.getSize()));
}
