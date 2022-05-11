#ifdef _USE_OMP_
#  include<omp.h>
#endif
#include<iostream>

/**
 * valgrind found some memory leaks using main-u1.cc
 * This is to check if the leaks also appear if something trivial is done with OpenMP
 * -> they do
 * */

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

//~ valgrind --leak-check=full --show-leak-kinds=all ./../try
//~ ==203229== Memcheck, a memory error detector
//~ ==203229== Copyright (C) 2002-2017, and GNU GPL'd, by Julian Seward et al.
//~ ==203229== Using Valgrind-3.15.0 and LibVEX; rerun with -h for copyright info
//~ ==203229== Command: ./../try
//~ ==203229== 
//~ 1.5
//~ 1.5
//~ ==203229== 
//~ ==203229== HEAP SUMMARY:
//~ ==203229==     in use at exit: 2,304 bytes in 5 blocks
//~ ==203229==   total heap usage: 9 allocs, 4 frees, 108,864 bytes allocated
//~ ==203229== 
//~ ==203229== 8 bytes in 1 blocks are still reachable in loss record 1 of 5
//~ ==203229==    at 0x483B7F3: malloc (in /usr/lib/x86_64-linux-gnu/valgrind/vgpreload_memcheck-amd64-linux.so)
//~ ==203229==    by 0x4A4E24C: ??? (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x4A5EBAA: ??? (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x4A4C679: ??? (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x4011B89: call_init.part.0 (dl-init.c:72)
//~ ==203229==    by 0x4011C90: call_init (dl-init.c:30)
//~ ==203229==    by 0x4011C90: _dl_init (dl-init.c:119)
//~ ==203229==    by 0x4001139: ??? (in /usr/lib/x86_64-linux-gnu/ld-2.31.so)
//~ ==203229== 
//~ ==203229== 24 bytes in 1 blocks are still reachable in loss record 2 of 5
//~ ==203229==    at 0x483B723: malloc (in /usr/lib/x86_64-linux-gnu/valgrind/vgpreload_memcheck-amd64-linux.so)
//~ ==203229==    by 0x483E017: realloc (in /usr/lib/x86_64-linux-gnu/valgrind/vgpreload_memcheck-amd64-linux.so)
//~ ==203229==    by 0x4A4E2AC: ??? (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x4A5D581: ??? (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x4A548E0: GOMP_parallel (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x10939F: tryreduce(int) (in /home/christiane/Documents/SU2-Carsten/su2/build/debug/try)
//~ ==203229==    by 0x109553: main (in /home/christiane/Documents/SU2-Carsten/su2/build/debug/try)
//~ ==203229== 
//~ ==203229== 192 bytes in 1 blocks are still reachable in loss record 3 of 5
//~ ==203229==    at 0x483B7F3: malloc (in /usr/lib/x86_64-linux-gnu/valgrind/vgpreload_memcheck-amd64-linux.so)
//~ ==203229==    by 0x4A4E24C: ??? (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x4A5C9A0: ??? (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x4A548C9: GOMP_parallel (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x10939F: tryreduce(int) (in /home/christiane/Documents/SU2-Carsten/su2/build/debug/try)
//~ ==203229==    by 0x109553: main (in /home/christiane/Documents/SU2-Carsten/su2/build/debug/try)
//~ ==203229== 
//~ ==203229== 304 bytes in 1 blocks are possibly lost in loss record 4 of 5
//~ ==203229==    at 0x483DD99: calloc (in /usr/lib/x86_64-linux-gnu/valgrind/vgpreload_memcheck-amd64-linux.so)
//~ ==203229==    by 0x40149CA: allocate_dtv (dl-tls.c:286)
//~ ==203229==    by 0x40149CA: _dl_allocate_tls (dl-tls.c:532)
//~ ==203229==    by 0x4DF2322: allocate_stack (allocatestack.c:622)
//~ ==203229==    by 0x4DF2322: pthread_create@@GLIBC_2.2.5 (pthread_create.c:660)
//~ ==203229==    by 0x4A5CDEA: ??? (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x4A548E0: GOMP_parallel (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x10939F: tryreduce(int) (in /home/christiane/Documents/SU2-Carsten/su2/build/debug/try)
//~ ==203229==    by 0x109553: main (in /home/christiane/Documents/SU2-Carsten/su2/build/debug/try)
//~ ==203229== 
//~ ==203229== 1,776 bytes in 1 blocks are still reachable in loss record 5 of 5
//~ ==203229==    at 0x483B7F3: malloc (in /usr/lib/x86_64-linux-gnu/valgrind/vgpreload_memcheck-amd64-linux.so)
//~ ==203229==    by 0x4A4E24C: ??? (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x4A5C7FB: ??? (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x4A548C9: GOMP_parallel (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
//~ ==203229==    by 0x10939F: tryreduce(int) (in /home/christiane/Documents/SU2-Carsten/su2/build/debug/try)
//~ ==203229==    by 0x109553: main (in /home/christiane/Documents/SU2-Carsten/su2/build/debug/try)
//~ ==203229== 
//~ ==203229== LEAK SUMMARY:
//~ ==203229==    definitely lost: 0 bytes in 0 blocks
//~ ==203229==    indirectly lost: 0 bytes in 0 blocks
//~ ==203229==      possibly lost: 304 bytes in 1 blocks
//~ ==203229==    still reachable: 2,000 bytes in 4 blocks
//~ ==203229==         suppressed: 0 bytes in 0 blocks
//~ ==203229== 
//~ ==203229== For lists of detected and suppressed errors, rerun with: -s
//~ ==203229== ERROR SUMMARY: 1 errors from 1 contexts (suppressed: 0 from 0)
