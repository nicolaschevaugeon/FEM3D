#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include "omp.h"



int main(){
#pragma omp parallel
  {
    //int num = 0;
    //#pragma omp for
    //for(size_t i = 0; i < 64; ++i) num += rand();
    //auto cpu = sched_getcpu();
    std::ostringstream os;
    os<<"Thread "<<omp_get_thread_num()<<" on cpu "<<sched_getcpu()<<std::endl;
    std::cout<<os.str()<<std::flush;
    
  }

  {
    size_t m = 1<<16;
    
    int nloop = 1<<16;
    std::vector<double > a(m);
    std::vector<double > b(m);
    std::vector<double > c(m);
    //auto now = std::chrono::high_resolution_clock::now;
    //typedef std::chrono::duration<double, std::milli> dt_t;
    //auto tstart = now();
    auto tstart= omp_get_wtime();
    for(int k = 0; k < nloop; ++k)
    for(size_t i= 0; i < m; ++i ){
      a[i] = cos(a[i]) + sin(b[i]);
      a[i] = log(abs(c[i]));
      a[i] = cos( sin(c[i] + sin(a[i]) + sqrt(fabs(a[i]))));
      c[i] = cos(a[i]) + sin(b[i]);
      c[i] = log(abs(c[i]));
      c[i] = cos( sin(c[i] + sin(a[i]) + sqrt(fabs(a[i]))));
      a[i] = cos(a[i]) + sin(b[i]);
      a[i] = log(abs(c[i]));
      a[i] = cos( sin(c[i] + sin(a[i]) + sqrt(fabs(a[i]))));
      c[i] = cos(a[i]) + sin(b[i]);
      c[i] = log(abs(c[i]));
      c[i] = cos( sin(c[i] + sin(a[i]) + sqrt(fabs(a[i]))));
    }
    auto tend = omp_get_wtime();
    auto dt  = tend-tstart;
    //std::cout <<dt.count() << std::endl;
    std::cout << dt  << std::endl;

    // tstart = now();
    tstart= omp_get_wtime();
    
#pragma omp parallel    
    for(int k = 0; k < nloop; ++k){
      /* std::ostringstream os;
      os<<"Thread "<<omp_get_thread_num()<<" on cpu "<<sched_getcpu()<<std::endl;
      std::cout<<os.str()<<std::flush;*/
#pragma omp for schedule( static, 128) 
      for(size_t i= 0; i < m; ++i ){
	a[i] = cos(a[i]) + sin(b[i]);
	a[i] = log(abs(c[i]));
	a[i] = cos( sin(c[i] + sin(a[i]) + sqrt(fabs(a[i]))));
	c[i] = cos(a[i]) + sin(b[i]);
	c[i] = log(abs(c[i]));
	c[i] = cos( sin(c[i] + sin(a[i]) + sqrt(fabs(a[i]))));
	a[i] = cos(a[i]) + sin(b[i]);
	a[i] = log(abs(c[i]));
	a[i] = cos( sin(c[i] + sin(a[i]) + sqrt(fabs(a[i]))));
	c[i] = cos(a[i]) + sin(b[i]);
	c[i] = log(abs(c[i]));
	c[i] = cos( sin(c[i] + sin(a[i]) + sqrt(fabs(a[i]))));
      } 
    }
    //tend = now();
    tend = omp_get_wtime();
    auto  dtpar  = tend-tstart;
    std::cout << dtpar << std::endl;
    std::cout << dt/dtpar << std::endl;
  }
}  
