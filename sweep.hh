#pragma once

#include"gaugeconfig.hh"
#include"gaugeconfig.hh"
#include"accum_type.hh"
#include"random_element.hh"
#include"get_staples.hh"

#ifdef _USE_OMP_
#  include<omp.h>
#endif

#include<random>
#include<vector>

template<class URNG, class Group> double sweep(gaugeconfig<Group> &U, vector<URNG> engine,
                                               const double delta, 
                                               const size_t N_hit, const double beta, const double xi=1.0, bool anisotropic=false) {

  std::uniform_real_distribution<double> uniform(0., 1.);
  typedef typename accum_type<Group>::type accum;
  size_t rate = 0;
  #ifdef _USE_OMP_
  int threads = omp_get_max_threads();
  double * omp_acc = new double[threads];
  #pragma omp parallel
  {
    int thread_num = omp_get_thread_num();
  #else
  int thread_num=0;  
  #endif  
  size_t temp=0;
  #pragma omp for
  for(size_t x0 = 0; x0 < U.getLt(); x0+=2) {
//OpenMP does not allow loop declaration as it was done for the other dimensions, still have to figure out why
    Group R;
    for(size_t x1 = 0; x1 < U.getLx(); x1++) {
      for(size_t x2 = 0; x2 < U.getLy(); x2++) {
        for(size_t x3 = 0; x3 < U.getLz(); x3++) {
            std::vector<size_t> x = {x0, x1, x2, x3};
          for(size_t mu = 0; mu < U.getndims(); mu++) {
            accum K;
            get_staples(K, U, x, mu, xi, anisotropic);
            for(size_t n = 0; n < N_hit; n++) {
              random_element(R, engine[thread_num], delta);
              double deltaS = beta/static_cast<double>(U.getNc())*
                (retrace(U(x, mu) * K) - retrace(U(x, mu) * R * K));
              bool accept = (deltaS < 0);
              if(!accept) accept = (uniform(engine[thread_num]) < exp(-deltaS));
              if(accept) {
                U(x, mu) = U(x, mu) * R;
                U(x, mu).restoreSU();
                temp += 1;
              }
            }
          }
        }
      }
    }
  }
  #pragma omp for
  for(size_t x0 = 1; x0 < U.getLt(); x0+=2) {
//OpenMP does not allow loop declaration as it was done for the other dimensions, still have to figure out why
    Group R;
    for(size_t x1 = 0; x1 < U.getLx(); x1++) {
      for(size_t x2 = 0; x2 < U.getLy(); x2++) {
        for(size_t x3 = 0; x3 < U.getLz(); x3++) {
            std::vector<size_t> x = {x0, x1, x2, x3};
          for(size_t mu = 0; mu < U.getndims(); mu++) {
            accum K;
            get_staples(K, U, x, mu, xi, anisotropic);
            for(size_t n = 0; n < N_hit; n++) {
              random_element(R, engine[thread_num], delta);
              double deltaS = beta/static_cast<double>(U.getNc())*
                (retrace(U(x, mu) * K) - retrace(U(x, mu) * R * K));
              bool accept = (deltaS < 0);
              if(!accept) accept = (uniform(engine[thread_num]) < exp(-deltaS));
              if(accept) {
                U(x, mu) = U(x, mu) * R;
                U(x, mu).restoreSU();
                temp += 1;
              }
            }
          }
        }
      }
    }
  }
  //maybe try out reduction(+: temp) in pragma at some point? But still need if-construct to close bracket
  #ifdef _USE_OMP_
    omp_acc[thread_num] = temp;
    rate = 0.;
  }
  for(size_t i = 0; i < threads; i++) {
    rate += omp_acc[i];
  }
  delete [] omp_acc;
  #else
  rate = temp;
  #endif
  return( double(rate)/double(N_hit)/double(U.getSize()));
}

