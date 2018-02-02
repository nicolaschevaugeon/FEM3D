#ifndef _fem_test_
#define _fem_test_

#include "fem.h"



template <class Scal, class Law >
  void oneelemtest(const coordinate<Scal > &X0, const coordinate<Scal > &X1,    const coordinate<Scal > &X2,  const coordinate<Scal > &X3, const coordinate<Scal > &x0, const coordinate<Scal > &x1,    const coordinate<Scal > &x2,  const coordinate<Scal > &x3, const Law &law ){
 
  using scal_t = Scal;
  using fem_t = fem< Scal >;
  fem_t pb;
  classification classif = {3, 10};
  pb.coordinates[0].resize(4);
  pb.coordinates[1].resize(4);
  pb.coordinates[0][0] = X0;
  pb.coordinates[0][1] = X1;
  pb.coordinates[0][2] = X2;
  pb.coordinates[0][3] = X3;

  pb.coordinates[1][0] = x0;
  pb.coordinates[1][1] = x1;
  pb.coordinates[1][2] = x2;
  pb.coordinates[1][3] = x3;

  pb.vertices ={ { &classif , 1}, { &classif , 2}, { &classif , 3 },  { &classif , 0 }  };
  pb.tets = { {&classif, { &pb.vertices[0],  &pb.vertices[1],  &pb.vertices[2], &pb.vertices[3] } , 0 }  };
  
  auto F = pb.F(pb.tets[0], 0, 1);
 
  auto P = law.PKI(F);
  scal_t psi = law.elastic_potential(F);
  std::cout << "F " << F << std::endl;
  std::cout << "F " << pb.Fv2(pb.tets[0], 0, 1)  << std::endl;
  std::cout << "P " << P << std::endl;
  std::cout << "psi " << P << std::endl;
  law_test(F, law);
 
  mat34cm< scal_t > _fint;
  auto & fint = _fint.data;
  scal_t eps = 0.00001; 
  for (size_t a=0; a< 4; ++a){
    for (size_t i=0; i< 3; ++i){
      int node = pb.tets[0].vertices[a]->index;
      pb.coordinates[1][node].data[i]+=eps;
      auto   Feps = pb.F(pb.tets[0], 0, 1);
      std::cout << i << " dF " << (Feps-F)/eps<< std::endl;
      std::cout << "Fint"<<i<<a << " " << contract2((Feps-F)/eps,P) <<std::endl;
      scal_t psieps = law.elastic_potential(Feps);
      fint[i+3*a] = (psieps-psi)/eps;
      // std::cout << psieps<< " " << psi << std::endl; 
      pb.coordinates[1][node].data[i]-=eps;
    }
  }
  
  std::cout << " fint num " << _fint << std::endl;

 
  auto op =  pb.getFop(pb.tets[0], 0  );
  std::cout << " op " << op  << std::endl; 
  _fint  = apply(op, P);
  std::cout << " fint op  " << _fint << std::endl;
  mat34cm<Scal >dx00 = {0.,-1.,0., 0.,0.,0., 0.,0.,0., 0.,0.,0.};
  
  
 
    

  mat34cm<Scal > fint2 = {0.,0.,0., 0.,0.,0., 0.,0.,0., 0.,0.,0.};
  tensor2_9cm<Scal > basis0 =  {pb.get_basis(pb.tets[0], 0)};
  for (int i = 0; i < 3; ++i){
    tensor2_9cm<Scal > dbasis1 = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
    dbasis1.data[i] = -1.;
    dbasis1.data[3+i] = -1.;
    dbasis1.data[6+i] = -1.;
    fint2.data[i] =  contract2(dbasis1*invert(basis0), P);
    std::cout << "i "<< " dF " <<  dbasis1*invert(basis0) << std::endl;
  }
  for (int i = 3; i < 12; ++i){
    tensor2_9cm<Scal > dbasis1 = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
    dbasis1.data[i-3] = 1.;
    fint2.data[i] =  contract2(dbasis1*invert(basis0), P);
    std::cout << "i "<< " dF " <<  dbasis1*invert(basis0) << std::endl;
  }
  std::cout << " fint slow  " <<fint2 << std::endl;
  



}



#endif
