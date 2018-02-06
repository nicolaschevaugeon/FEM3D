#include "tensor2_9cm.h"
#include "hyperelastic.h"

int main(){
   
  using scal_t = double;
  tensor2_9cm< scal_t> F = {-4.,2.,3.,1,4.,2.,1.,5.,2.};

  Hyperelastic_StVenantKirchhoff < scal_t > lawst = {1.d,2.d};
  Hyperelastic_NeoHookeen < scal_t > lawneo = {1.d,2.d};
  
  std::cout << "testing lawneo" << std::endl;
  bool failed = law_test(F, lawneo);
  if (failed){
    std::cout  << "Test on law  Hyperelastic_NeoHookeen Failed " << std::endl;
    return 1;
  }
  std::cout << "testing lawstvenant" << std::endl;
  failed = law_test(F, lawst);
  if (failed){
    std::cout  << "Test on law   Hyperelastic_StVenantKirchhoff Failed " << std::endl;
    return 1;
  }
  return 0;
}

