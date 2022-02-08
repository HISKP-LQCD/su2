#ifdef _USE_OMP_
#  include<omp.h>
#endif
#include<iostream>

double tryreduce(const int a){
  double sum=0;
  //try to add up numbers by reduction operation
  #pragma omp parallel for reduction(+:sum)
  for(int i=0;i<a;i+=1){
    sum+=0.1*i;
  }
  std::cout << sum << std::endl;
  
  sum=0;
  //try to add up numbers through array accessed with thread number
  #ifdef _USE_OMP_
  int threads = omp_get_max_threads();
  static double * omp_acc = new double[threads];
  #pragma omp parallel
  {
    int thread_num = omp_get_thread_num();
  #endif
  double tmp=0;
  #pragma omp for
  for(int i=0; i<a; i+=1){
    tmp+=0.1*i;
  }
  #ifdef _USE_OMP_
    omp_acc[thread_num] = tmp;
  }
  for(size_t i = 0; i < threads; i++) {
    sum += omp_acc[i];
  }
  delete [] omp_acc;
  #else
  sum += tmp;
  #endif
  std::cout << sum << std::endl;
  return sum;
}


int main(int ac, char* av[]) {
    tryreduce(6);
    return(0);
}
