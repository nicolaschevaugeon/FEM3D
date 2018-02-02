#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>

struct tensor{
  double data[9];
};

tensor operator+(const tensor & a, const tensor &b){
  return tensor{ a.data[0]+b.data[0], a.data[1]+b.data[1], a.data[2]+b.data[2],a.data[3]+b.data[3],a.data[4]+b.data[4],a.data[5]+b.data[5],a.data[6]+b.data[6], a.data[7]+b.data[7], a.data[8]+b.data[8]   };
}

struct tensor_v{
  typedef double  d16si __attribute__((vector_size(4*sizeof(double) )));
  d16si data[3];
};

tensor_v operator+(const tensor_v & a, const tensor_v &b){
  tensor_v r  = {a.data[0] +b.data[0], a.data[1]+b.data[1], a.data[2]+b.data[2] };
  return r;
}


int main()
{
  int m = 1<< 24;
  {
    std::vector <tensor> A(m);
    std::vector <tensor> B(m);
    std::vector <tensor> C(m);
    //auto t_start = std::chrono::high_resolution_clock::now();
    std::transform( A.begin(), A.end(), B.begin(), C.begin(), [](const tensor & ai, const tensor & bi){return ai+bi;});
    //auto t_end = std::chrono::high_resolution_clock::now();
    //std::cout <<  "Tensor Add  transform " << " " <<  std::chrono::duration<double, std::milli>(t_end-t_start).count() << std::endl;
  }
  /* {
    std::vector <tensor_v> A(m);
    std::vector <tensor_v> B(m);
    std::vector <tensor_v> C(m);
    auto t_start = std::chrono::high_resolution_clock::now();
    std::transform( A.begin(), A.end(), B.begin(), C.begin(), [](const tensor_v & ai, const tensor_v & bi){return ai+bi;});
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout <<  "Tensor Add  transform " << " " <<  std::chrono::duration<double, std::milli>(t_end-t_start).count() << std::endl;
    }*/
  
  
  
}
     
