#ifndef  _hyperelastic_
#define  _hyperelastic_
#include "tensor2_9cm.h"
#include <iostream>

template < class Scal >
class Hyperelastic{
public:
  virtual tensor2_9cm< Scal > PKI(const tensor2_9cm< Scal > & F ) const = 0;
};


template <class Scal >
class Hyperelastic_StVenantKirchhoff : public Hyperelastic< Scal >{
public:
  Hyperelastic_StVenantKirchhoff(Scal _lambda, Scal _mu):lambda(_lambda), mu(_mu){}
  tensor2_9cm< Scal > PKI(const tensor2_9cm< Scal > & F ) const{
    return PKI<tensor2_9cm< Scal > > (F);
  }
  template < class Tensor2 >
  Tensor2 PKI(const Tensor2 & F ) const {
    const auto E = Green_Lagrange_sym(F);
    return  F*( tensor2_isotrope< typename Tensor2::value_type >{ lambda*trace(E) } +(2.f*mu)*E);
  }
  template < class Tensor2 >
  Scal elastic_potential(const Tensor2 & F ) const {
    const auto E = Green_Lagrange_sym(F);
    const Scal trE =  trace(E);
    return 0.5*lambda*trE*trE + mu*I2(E);
  }
  /*
  template < class Tensor2 >
  Scal elastic_potentialE(const Tensor2 & E ) const {
    const Scal trE = trace(E);
    const Scal trE2 = I2(E);
    return 0.5*lambda*trE*trE + mu*trE2;
    }*/
  const Scal lambda, mu;
};

template <class Scal >
class Hyperelastic_NeoHookeen : public Hyperelastic < Scal >{
public:
  Hyperelastic_NeoHookeen(Scal _lambda, Scal _mu):lambda(_lambda), mu(_mu){}
  tensor2_9cm< Scal > PKI(const tensor2_9cm< Scal > & F ) const{
    return PKI<tensor2_9cm< Scal > > (F);
  }
  template < class Tensor2 >
  Tensor2 PKI(const Tensor2 & F ) const {
    const Scal lnJ = log(det(F));
    const auto invC = invert(Right_Cauchy_Green_sym(F));
    return F*( tensor2_isotrope< typename Tensor2::value_type >{ mu } +  (lambda*lnJ-mu)*invC );
    
  }
  template < class Tensor2 >
  Scal elastic_potential(const Tensor2 & F ) const {
    const Scal lnJ = log(det(F));
    const auto C = Right_Cauchy_Green_sym(F);
    const Scal trC = trace(C);
    return mu*(0.5*trC-1.5-lnJ) +  lambda*0.5*lnJ*lnJ;
  }
  const Scal lambda, mu;
};


template <class Scal, class Law > void law_test( const tensor2_9cm<Scal > &F, const Law & law ){
  const tensor2_9cm<Scal > P = law.PKI(F);
  const Scal psi = law.elastic_potential(F);
  //std::cout << "psi "<< psi << std::endl;
  // throw;
  tensor2_9cm<Scal > numP;
  const Scal eps = 0.01;
  for (int i = 0; i< 9; ++i){
    tensor2_9cm<Scal > Fpeps = F;
    tensor2_9cm<Scal > Fmeps = F;
    Fpeps.data[i] += eps;
    Fmeps.data[i] -= eps;
    const Scal Pi = (law.elastic_potential(Fpeps) - law.elastic_potential(Fmeps))/2./eps;
    numP.data[i] = Pi;
  }
  if( norm2(P-numP) > eps*norm2(numP)){
    std::cout << "Error in law_test " << __FILE__<< ":"<< __LINE__<< std::endl;
    std::cout << "P    " << P   << std::endl;
    std::cout << "dfP  " << numP << std::endl;
    std::cout << norm2(P-numP)*norm2(numP) << " "<<  eps*norm2(numP) << std::endl;
    throw;
  }
}

#endif
