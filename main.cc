#include <array>
#include <vector>
#include <parallel/algorithm>
#include <chrono>
#include <iostream>
#include <omp.h>
#include "tensor2_9cm.h"
#include <memory>
#include<functional>
#include<fstream>
#include <regex>
#include "gmsh_import.h"
#include "hyperelastic.h"
#include "fem.h"
#include "fem_test.h"


template < class ITER1, class ITER2, class OP>
void mytransform(const ITER1 &b1, const ITER1 &e1 , const ITER2 &b2, const OP& op){
  size_t nth = omp_get_max_threads();
  std::vector < ITER1> b1_t(nth);
  std::vector < ITER1> e1_t(nth);
  std::vector < ITER2> b2_t(nth);
  const size_t l = std::distance(b1,e1);
  for (size_t i = 0; i < nth; ++i){
    b1_t[i] = b1+i*l/nth;
    b2_t[i] = b2+i*l/nth;
    e1_t[i] = b1+(i+1)*l/nth;
  }
  e1_t[nth-1] = e1;
#pragma omp parallel for
  for (size_t i =0; i< nth; ++i)
    std::transform(b1_t[i], e1_t[i],  b2_t[i], op  );
}



int main(int argc, char *argv[]){
  if (argc < 2){
    std::cerr  << "Usage : " << argv[0] << " meshfilename.msh [numthreads]"<< std::endl;
    return 1;
  }

  
  int numthreads = 4 ;
  if (argc == 3) numthreads = std::stoi (argv[2]);
  omp_set_num_threads(numthreads);
  
  using scal_t = double;
  tensor2_9cm< scal_t> F = {-4.,2.,3.,1,4.,2.,1.,5.,2.};

  Hyperelastic_StVenantKirchhoff < scal_t > lawst = {1.d,2.d};
  Hyperelastic_NeoHookeen < scal_t > lawneo = {1.d,2.d};
  
  std::cout << "testing lawneo" << std::endl;
  law_test(F, lawneo);
  std::cout << "testing lawstvenant" << std::endl;
  law_test(F, lawst);

 
  /*tensor2_9cm<double > A = {1.,2.,3.,1,4.,2.,1.,5.,2.};
  std::cout << invert(A) << std::endl;
  std::cout << invert(A)*A << std::endl;
  std::cout << A << std::endl;
  std::cout << transpose(A) << std::endl;
  std::cout << transpose(A)*transpose(invert(A)) << std::endl;
  std::cout << invert_transpose(A)*transpose(A) << std::endl;
  */

  //std::fstream inputfilegmsh("cube_20.msh");
  //std::fstream inputfilegmsh("cube_100.msh");
  std::string meshfilename = argv[1];
  std::ifstream inputfilegmsh(meshfilename);
  if(!inputfilegmsh.is_open() ){
    std::cerr << "Can't open meshfile " <<  meshfilename << std::endl;
    return 1;
  }
  gmsh_import::mesh m = gmsh_import::import_GMSH(inputfilegmsh);
  inputfilegmsh.close();
  std::cout << "Read gmsh mesh, " << "nv =" << m.nodes.size() << "nelem = " << m.elems.size() << std::endl;
  
  
  using fem_t = fem< scal_t >;
  fem_t pb;

  // Putting the mesh in fem
  pb.coordinates[0].resize(m.nodes.size());
  for (size_t i=0; i < m.nodes.size(); ++i){
    pb.coordinates[0][i] = {static_cast< scal_t> (m.nodes[i].x), static_cast< scal_t> ( m.nodes[i].y), static_cast< scal_t> (m.nodes[i].z)};
  }
  pb.coordinates[1] =  pb.coordinates[0];
  pb.vertices.resize(m.nodes.size() );
  classification classif = {3, 10};
  for( size_t i = 0; i < pb.vertices.size(); ++i){
    pb.vertices[i] = {&classif, i};
  }
  size_t ntets = 0;
  for (size_t i=0; i < m.elems.size(); ++i){
    if (m.elems[i].type == 4){
      auto tet= m.elems[i];
      int k0= tet.gmshnodeids[0];
      int k1= tet.gmshnodeids[1];
      int k2= tet.gmshnodeids[2];
      int k3= tet.gmshnodeids[3];
      pb.tets.push_back( fem_t::tetrahedre {&classif, { &pb.vertices[k0-1],  &pb.vertices[k1-1], &pb.vertices[k2-1], &pb.vertices[k3-1] }, ntets } );
      ++ntets;
    }
  }
  m.clear(); // I don't need gmsh mesh representation anymore ...
  
  // setting deformed nodes coordinates
   for_each( pb.coordinates[1].begin(), pb.coordinates[1].end(), []( coordinate< scal_t> &X  ){
       const coordinate<scal_t> u ={2.,3.,4.};
       const tensor2_9cm< scal_t> F = {1.3, 0.5, 0.2,1.2, 1.01, 0.6, 2., 0.7, 0.37 };
       X =   F*X +u;
    }

    );
  
  std::cout << "NTETS " << ntets;
  std::cout << pb.tets[0] << std::endl;
 
 
 

  const Hyperelastic< scal_t > & law = lawneo;
  pb.law = &law;
  
  int ne = ntets;
 
  __gnu_parallel::_Parallelism par_tag = __gnu_parallel::_Parallelism::parallel_unbalanced; //130

  pb.setFops();
 
  auto now = std::chrono::high_resolution_clock::now;
  typedef std::chrono::duration<double, std::milli> dt_t;
  dt_t time_serial;
  dt_t time_parallel;
  
  std::cout << "Computing F, v1 "<< std::endl;
  
  for (int i = 0; i < 1; ++i)
    {
      std::vector< fem_t::tensor_t > Fs(ne);
      auto t_start = now();
      __gnu_parallel::transform(pb.tets.begin(), pb.tets.end(), Fs.begin(),  [ &pb ](const fem_t::tetrahedre & tet ){ return pb.F(tet,0,1);  }, par_tag);
      auto t_end = now();
      time_parallel = t_end-t_start;
      //std::cout << Fs[0] << std::endl;
      //std::cout << Fs[ne-1] << std::endl;
      
     }
  for (int i = 0; i < 1; ++i)
    {
      std::vector< fem_t::tensor_t > Fs(ne);
      auto t_start = now();
      std::transform(pb.tets.begin(), pb.tets.end(), Fs.begin(),  [&pb](const fem_t::tetrahedre & tet ){ return pb.Fv2(tet,0,1);  });
      auto t_end = now();
      time_serial = t_end-t_start;
      //std::cout << Fs[0] << std::endl;
      //std::cout << Fs[ne-1] << std::endl;
    }
  std::cout <<  " parallel transform  "<< omp_get_max_threads() << " " << time_parallel.count() << std::endl;
  std::cout <<  " serial   transform  "<< 1  << " " << time_serial.count() << std::endl;
  std::cout <<  " speed up "<<  " " << time_serial.count()/time_parallel.count() << std::endl;

  std::cout << "Computing F, v2 "<< std::endl;
  
 for (int i = 0; i < 1; ++i)
    {
      std::vector< fem_t::tensor_t > Fs(ne);
      auto t_start = now();
      __gnu_parallel::transform(pb.tets.begin(), pb.tets.end(), Fs.begin(),  [ &pb ](const fem_t::tetrahedre & tet ){ return pb.Fv2(tet,0,1);  }, par_tag);
      auto t_end = now();
      time_parallel = t_end-t_start;
      // std::cout << Fs[0] << std::endl;
      //std::cout << Fs[ne-1] << std::endl;
      
     }

 for (int i = 0; i < 1; ++i)
   {
     std::vector< fem_t::tensor_t > Fs(ne);
     auto t_start = now();
     std::transform(pb.tets.begin(), pb.tets.end(), Fs.begin(),  [&pb](const fem_t::tetrahedre & tet ){ return pb.F(tet,0,1);  });
      auto t_end = now();
       time_serial = t_end-t_start;
   }
        

 std::cout <<  " parallel transform  "<< omp_get_max_threads() << " " << time_parallel.count() << std::endl;
  std::cout <<  " serial   transform  "<< 1  << " " << time_serial.count() << std::endl;
  std::cout <<  " speed up "<<  " " << time_serial.count()/time_parallel.count() << std::endl;

  std::cout << "Computing F, v3 (using op) "<< std::endl;
 for (int i = 0; i < 1; ++i)
    {
      std::vector< fem_t::tensor_t > Fs(ne);
      auto t_start = now();
      __gnu_parallel::transform(pb.tets.begin(), pb.tets.end(), Fs.begin(),  [ &pb ](const fem_t::tetrahedre & tet ){ return pb.Fv3(tet,1);  }, par_tag);
      auto t_end = now();
      time_parallel = t_end-t_start;
    }

 for (int i = 0; i < 1; ++i)
   {
     std::vector< fem_t::tensor_t > Fs(ne);
     auto t_start = now();
     std::transform(pb.tets.begin(), pb.tets.end(), Fs.begin(),  [&pb](const fem_t::tetrahedre & tet ){ return pb.Fv3(tet,1);  });
     auto t_end = now();
     time_serial = t_end-t_start;
   }
 std::cout <<  " parallel transform  "<< omp_get_max_threads() << " " << time_parallel.count() << std::endl;
 std::cout <<  " serial   transform  "<< 1  << " " << time_serial.count() << std::endl;
 std::cout <<  " speed up "<<  " " << time_serial.count()/time_parallel.count() << std::endl;
 
 /// 
  std::cout << "Computing " << ne << " P (F on the fly, law neo)" << std::endl;
  for (int i = 0; i < 1; ++i)
       {
	 std::vector< fem_t::tensor_t >  Ps(ne);
	 auto t_start = now();
	 __gnu_parallel::transform(pb.tets.begin(), pb.tets.end(), Ps.begin(),  [&pb, &lawneo](const fem_t::tetrahedre & tet ){ return lawneo.PKI(pb.F(tet,0,1));  }, par_tag);
	 
	 auto t_end = now();
	 time_parallel = t_end-t_start;
       }
  
  for (int i = 0; i < 1; ++i)
    {
      std::vector< fem_t::tensor_t > Ps(ne);
      auto t_start = now();
      std::transform(pb.tets.begin(), pb.tets.end(), Ps.begin(),  [&pb, &lawneo](const fem_t::tetrahedre & tet ){ return lawneo.PKI(pb.F(tet,0,1));  });
      auto t_end = now();
      time_serial = t_end-t_start;
    }
  std::cout <<  " parallel transform  "<< omp_get_max_threads() << " " << time_parallel.count() << std::endl;
  std::cout <<  " serial   transform  "<< 1  << " " << time_serial.count() << std::endl;
  std::cout <<  " speed up "<<  " " << time_serial.count()/time_parallel.count() << std::endl;

  std::cout << "Computing " << ne << " P (F on the fly, base law )" << std::endl;
  for (int i = 0; i < 1;  ++i)
       {
	 std::vector< fem_t::tensor_t >  Ps(ne);
	 auto t_start = now();
	 __gnu_parallel::transform(pb.tets.begin(), pb.tets.end(), Ps.begin(),  [&pb, &law](const fem_t::tetrahedre & tet ){ return law.PKI(pb.F(tet,0,1));  }, par_tag);
	 
	 auto t_end = now();
	 time_parallel = t_end-t_start;
       }
  
  for (int i = 0; i < 1; ++i)
    {
      std::vector< fem_t::tensor_t > Ps(ne);
      auto t_start = now();
      std::transform(pb.tets.begin(), pb.tets.end(), Ps.begin(),  [&pb, &law](const fem_t::tetrahedre & tet ){ return law.PKI(pb.F(tet,0,1));  });
      auto t_end = now();
      time_serial = t_end-t_start;
    }
  std::cout <<  " parallel transform  "<< omp_get_max_threads() << " " << time_parallel.count() << std::endl;
  std::cout <<  " serial   transform  "<< 1  << " " << time_serial.count() << std::endl;
  std::cout <<  " speed up "<<  " " << time_serial.count()/time_parallel.count() << std::endl;


  // computing P (Fops)
  std::cout << "Computing " << ne << " P ( using precomputed op, law neo )" << std::endl;
  for (int i = 0; i < 1; ++i)
       {
	 std::vector< fem_t::tensor_t >  Ps(ne);
	 auto t_start = now();
	 __gnu_parallel::transform(pb.tets.begin(), pb.tets.end(), Ps.begin(),  [&pb, &lawneo](const fem_t::tetrahedre & tet ){ return lawneo.PKI(pb.Fv3(tet,1));  }, par_tag);
	 auto t_end = now();
	 time_parallel = t_end-t_start;
       }
  
  for (int i = 0; i < 1; ++i)
    {
      std::vector< fem_t::tensor_t > Ps(ne);
      auto t_start = now();
      std::transform(pb.tets.begin(), pb.tets.end(), Ps.begin(),  [&pb, &lawneo](const fem_t::tetrahedre & tet ){ return lawneo.PKI(pb.Fv3(tet,1));  });
      auto t_end = now();
      time_serial = t_end-t_start;
    }
  
  std::cout <<  " parallel transform  "<< omp_get_max_threads() << " " << time_parallel.count() << std::endl;
  std::cout <<  " serial   transform  "<< 1  << " " << time_serial.count() << std::endl;
  std::cout <<  " speed up "<<  " " << time_serial.count()/time_parallel.count() << std::endl;

  std::cout << "Computing " << ne << " P ( using precomputed op, base law )" << std::endl;
  for (int i = 0; i < 1;  ++i)
       {
	 std::vector< fem_t::tensor_t >  Ps(ne);
	 auto t_start = now();
	 __gnu_parallel::transform(pb.tets.begin(), pb.tets.end(), Ps.begin(),  [&pb, &law](const fem_t::tetrahedre & tet ){ return law.PKI(pb.Fv3(tet,1));  }, par_tag);
	 auto t_end = now();
	 time_parallel = t_end-t_start;
       }
  
  for (int i = 0; i < 1; ++i)
    {
      std::vector< fem_t::tensor_t > Ps(ne);
      auto t_start = now();
      std::transform(pb.tets.begin(), pb.tets.end(), Ps.begin(),  [&pb, &law](const fem_t::tetrahedre & tet ){ return law.PKI(pb.Fv3(tet,1));  });
      auto t_end = now();
      time_serial = t_end-t_start;
    }
  std::cout <<  " parallel transform  "<< omp_get_max_threads() << " " << time_parallel.count() << std::endl;
  std::cout <<  " serial   transform  "<< 1  << " " << time_serial.count() << std::endl;
  std::cout <<  " speed up "<<  " " << time_serial.count()/time_parallel.count() << std::endl;

  
  // computing Fint
  std::cout << "Computing Fint nelem = " << ne << " ( using precomputed op and law neo)" << std::endl;
  for (int i = 0; i < 1; ++i)
       {
	 std::vector< mat34cm< scal_t > >  Fint(ne);
	 auto t_start = now();
	 __gnu_parallel::transform(pb.tets.begin(), pb.tets.end(), Fint.begin(),  [&pb, &lawneo](const fem_t::tetrahedre & tet ){
	     const auto P = lawneo.PKI(pb.Fv3(tet,1));
	     return apply( pb.Fops[tet.id], P);
	   },
	   par_tag);
	 
	 auto t_end = now();
	 time_parallel = t_end-t_start;
	 //std::cout << Fint[0] << std::endl;
	 //std::cout << Fint[ne-1] << std::endl;
       }
  
  for (int i = 0; i < 1; ++i)
    {
      std::vector< mat34cm< scal_t > >  Fint(ne);
      auto t_start = now();
      std::transform(pb.tets.begin(), pb.tets.end(), Fint.begin(),  [&pb, &lawneo](const fem_t::tetrahedre & tet ){
	  const auto P = lawneo.PKI(pb.Fv3(tet,1));
	  return apply( pb.Fops[tet.id], P);
	}
	);
      
      auto t_end = now();
      time_serial = t_end-t_start;
      //std::cout <<  "serial transform (lawneo )"<< omp_get_max_threads() << " " <<  std::chrono::duration<double, std::milli>(t_end-t_start).count() << std::endl;
      //std::cout << Fint[0] << std::endl;
      //std::cout << Fint[ne-1] << std::endl;
    }
  std::cout <<  " parallel transform ( lawneo ) "<< omp_get_max_threads() << " " << time_parallel.count() << std::endl;
  std::cout <<  " serial   transform ( lawneo ) "<< 1  << " " << time_serial.count() << std::endl;
  std::cout <<  " serial   speed up "<<  " " << time_serial.count()/time_parallel.count() << std::endl;

 

   std::cout << "Computing Fint nelem = " << ne << " ( using precomputed op and base law)" << std::endl;  
  
  for (int i = 0; i < 1; ++i)
       {
	 std::vector< mat34cm< scal_t > >  Fint(ne);
	 auto t_start = now();
	 __gnu_parallel::transform(pb.tets.begin(), pb.tets.end(), Fint.begin(),  [&pb, &law](const fem_t::tetrahedre & tet ){
	     const auto P = law.PKI(pb.Fv3(tet,1));
	     return apply( pb.Fops[tet.id], P);
	   },
	   par_tag);
	 auto t_end = now();
	 time_parallel = t_end-t_start;
	 //std::cout << Fint[0] << std::endl;
	 //std::cout << Fint[ne-1] << std::endl;
       }
  
  for (int i = 0; i < 1; ++i)
    {
      std::vector< mat34cm< scal_t > >  Fint(ne);
      auto t_start = now();
      std::transform(pb.tets.begin(), pb.tets.end(), Fint.begin(),  [&pb, &law](const fem_t::tetrahedre & tet ){
	  const auto P = law.PKI(pb.Fv3(tet,1));
	  return apply( pb.Fops[tet.id], P);
	}
	);
      
      auto t_end = now();
      time_serial = t_end-t_start;
      std::cout << Fint[0] << std::endl;
      std::cout << Fint[ne-1] << std::endl;
    }
  std::cout <<  " parallel transform (base law ) "<< omp_get_max_threads() << " " << time_parallel.count() << std::endl;
  std::cout <<  " serial   transform (base law ) "<< 1  << " " << time_serial.count() << std::endl;
  std::cout <<  " serial   speed up "<<  " " << time_serial.count()/time_parallel.count() << std::endl;
    std::cout << "Computing Fint nelem = " << ne << " ( using precomputed op and base law)" << std::endl;  

 std::cout << "Computing Fint nelem = " << ne << " ( using fem member and fem->law)" << std::endl;  
    
  for (int i = 0; i < 1; ++i)
       {
	 std::vector< mat34cm< scal_t > >  Fint(ne);
	 auto t_start = now();
	 __gnu_parallel::transform(pb.tets.begin(), pb.tets.end(), Fint.begin(),  [&pb](const fem_t::tetrahedre & tet ){
	     return pb.internal_forces( tet, 1);
	   },
	   par_tag);
	 auto t_end = now();
	 time_parallel = t_end-t_start;
	 std::cout << Fint[0] << std::endl;
	 std::cout << Fint[ne-1] << std::endl;
       }
  
  for (int i = 0; i < 1; ++i)
    {
      std::vector< mat34cm< scal_t > >  Fint(ne);
      auto t_start = now();
      std::transform(pb.tets.begin(), pb.tets.end(), Fint.begin(),  [&pb](const fem_t::tetrahedre & tet ){
	  return pb.internal_forces( tet, 1);
	}
	);
      
      auto t_end = now();
      time_serial = t_end-t_start;
      std::cout << Fint[0] << std::endl;
      std::cout << Fint[ne-1] << std::endl;
    }
  auto time_parallel2 = time_serial;
  for (int i = 0; i < 1; ++i)
    {
      std::vector< mat34cm< scal_t > >  Fint(ne);
      auto t_start = now();
      mytransform(pb.tets.begin(), pb.tets.end(), Fint.begin(),  [&pb](const fem_t::tetrahedre & tet ){
	  return pb.internal_forces( tet, 1);
	}
	);
      
      auto t_end = now();
      time_parallel2 = t_end-t_start;
      std::cout << Fint[0] << std::endl;
      std::cout << Fint[ne-1] << std::endl;
    }
  std::cout <<  " parallel transform (base law ) "<< omp_get_max_threads() << " " << time_parallel.count() << std::endl;
  std::cout <<  " serial   transform (base law ) "<< 1  << " " << time_serial.count() << std::endl;
  std::cout <<  " parallel   mytransform (base law ) "<< 1  << " " << time_parallel2.count() << std::endl;
  std::cout <<  " serial   speed up "<<  " " << time_serial.count()/time_parallel.count() << std::endl;
  std::cout <<  " my       speed up "<<  " " << time_serial.count()/time_parallel2.count() << std::endl;
   
 
  return 0;
   
}
