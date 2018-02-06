#include "tensor2_9cm.h"

int main(){
  tensor2_9cm<double > A = {1.,2.,3.,1,4.,2.,1.,5.,2.};
  std::cout << invert(A) << std::endl;
  std::cout << invert(A)*A << std::endl;
  std::cout << A << std::endl;
  std::cout << transpose(A) << std::endl;
  std::cout << transpose(A)*transpose(invert(A)) << std::endl;
  std::cout << invert_transpose(A)*transpose(A) << std::endl;
}
