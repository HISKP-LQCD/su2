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

/**
 * go through the entire lattice
 * do N_hit Metropolis-Updates of every link (accept proposed change with probability min(1, exp(-deltaS)))
 * update only changes one link, not the surrounding links, but for the action, every plaquette including this link is needed
 * -> sum up the unchanged part in staples, only has to be calculated once for every link
 * If the lattice is anisotropic with anisotropy xi, the action weights the temporal (including links in direction 0) 
 * and spatial links differently, this is implemented in get_staples.hh.
 * For the update, the nearest neighbour links have to be constant, so parallelization is not trivial.
 * It is done by first updating all even time slices and then all odd time slices.
 * The overall acceptance rate, as well as the acceptance rate of the temporal links, are measured and returned.
 * An update means replacing one link U -> R*U, where R is a random element. 
 * The acceptance rate can be tuned with delta, which determines the possible regions from which R is drawn.
 * Is it more efficient to use two of every variable for measuring the rates, or would it be better to use vectors for everything?
 * When using vectors, a new reduction directive would have to be declared for the vectors
 * omp_priv has to be initialized with the = operator, and I could not find a way to initialize it to zero this way
 * pragma omp declare reduction (+ : std::vector<double> : omp_out += omp_in) initializer (omp_priv = ?)
 * So for the moment I would stick with the two separate variables
 * For the normalization of the temporal rate: There is only one temporal link for each lattice point, so the normalization is done with U.getVolume()
 * With this definition, for xi=1 the acceptance rates only differ in the third significant digit, so this is correct
 * */

template<class URNG, class Group> std::vector<double> sweep(gaugeconfig<Group> &U, std::vector<URNG> engine,
                                               const double &delta, 
                                               const size_t &N_hit, const double &beta,
                                               const double &xi=1.0, const bool& anisotropic=false) {

  std::uniform_real_distribution<double> uniform(0., 1.);
  typedef typename accum_type<Group>::type accum;
  size_t rate = 0, rate_time = 0;
  #ifdef _USE_OMP_
  #pragma omp parallel
  {
    size_t thread_num = omp_get_thread_num();
  #else
  size_t thread_num=0;  
  #endif  
  #pragma omp for reduction (+: rate, rate_time)
  for(size_t x0 = 0; x0 < U.getLt(); x0+=2) {
//Cannot use elements of a vector as iteration variables in for-loop with OpenMP, so use dummy variables
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
                rate += 1;
                if(mu == 0){
                  rate_time += 1;
                }
              }
            }
          }
        }
      }
    }
  }
  #pragma omp for reduction (+: rate, rate_time)
  for(size_t x0 = 1; x0 < U.getLt(); x0+=2) {
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
                rate += 1;
                if(mu == 0){
                  rate_time += 1;
                }
              }
            }
          }
        }
      }
    }
  }
  #ifdef _USE_OMP_
  }
  #endif
  std::vector<double> res = { double(rate)/double(N_hit)/double(U.getSize()) , double(rate_time)/double(N_hit)/double(U.getVolume()) };
  return res;
}

