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






int main(int argc, char *argv[]){
  using scal_t = double;
  tensor2_9cm< scal_t> F = {-4.,2.,3.,1,4.,2.,1.,5.,2.};

  Hyperelastic_StVenantKirchhoff < scal_t > lawst = {1.d,2.d};
  Hyperelastic_NeoHookeen < scal_t > lawneo = {1.d,2.d};
  
  std::cout << "testing lawneo" << std::endl;
  law_test(F, lawneo);
  std::cout << "testing lawstvenant" << std::endl;
  law_test(F, lawst);
  
  const coordinate<scal_t > X0 = {1.,2.,0.3};
  const coordinate<scal_t > X1 ={1./3.,0.1,0.2};
  const coordinate<scal_t > X2 = {0.7,1.,0.2};
  const coordinate<scal_t > X3 = {0.1,0.3,1.};
  const coordinate<scal_t > x0 ={0.1,0.25,0.17};
  const coordinate<scal_t > x1 ={3.,0.15,0.16};
  const coordinate<scal_t > x2 ={0.37,2.,0.29};
  const coordinate<scal_t > x3 = {0.42,0.11,2.};

  /*
  const coordinate<scal_t > X0 = {0., 0., 0.};
  const coordinate<scal_t > X1 = {1., 0., 0.};
  const coordinate<scal_t > X2 = {0., 1., 0.};
  const coordinate<scal_t > X3 = {0,  0., 1.};
  const coordinate<scal_t > x0 = {0., 0., 0.};
  const coordinate<scal_t > x1 = {2., 0., 0.};
  const coordinate<scal_t > x2 = {0., 3., 0.};
  const coordinate<scal_t > x3 = {0,  0., 4.};
  */
  
  
  //oneelemtest(X0, X1, X2, X3, x0, x1, x2, x3, lawneo);
  
  /*tensor2_9cm<double > A = {1.,2.,3.,1,4.,2.,1.,5.,2.};
  std::cout << invert(A) << std::endl;
  std::cout << invert(A)*A << std::endl;
  std::cout << A << std::endl;
  std::cout << transpose(A) << std::endl;
  std::cout << transpose(A)*transpose(invert(A)) << std::endl;
  std::cout << invert_transpose(A)*transpose(A) << std::endl;
  */

 
  //    std::fstream inputfilegmsh("cube_100.msh");
  std::fstream inputfilegmsh("cube_200.msh");
  gmsh_import::mesh m = gmsh_import::import_GMSH(inputfilegmsh);
  std::cout << "Read gmsh mesh, " << "nv =" << m.nodes.size() << "nelem = " << m.elems.size() << std::endl;
  
  
  using fem_t = fem< scal_t >;
  fem_t pb;
  pb.coordinates[0].resize(m.nodes.size());
  for (int i=0; i < m.nodes.size(); ++i){
    pb.coordinates[0][i] = {static_cast< scal_t> (m.nodes[i].x), static_cast< scal_t> ( m.nodes[i].y), static_cast< scal_t> (m.nodes[i].z)};
  }


  
  pb.coordinates[1] =  pb.coordinates[0];
  for_each( pb.coordinates[1].begin(), pb.coordinates[1].end(), []( coordinate< scal_t> &in  ){
      scal_t in0 = in.data[0];
      scal_t in1 = in.data[1];
      scal_t in2 = in.data[2];
      scal_t u0  =2.;
      scal_t u1  =3.;
      scal_t u2  =4.;
      scal_t F00 = 1.3;
      scal_t F01 = 1.2;
      scal_t F02 = 2.;
      scal_t F10 = 0.5;
      scal_t F11 = 1.01;
      scal_t F12 = 0.7;
      scal_t F20 = 0.2;
      scal_t F21 = 0.6;
      scal_t F22 = 0.37;
      scal_t out0 = F00*in0 +F01*in1+ F02*in2 + u0;
      scal_t out1 = F10*in0 +F11*in1+ F12*in2 + u1;
      scal_t out2 = F20*in0 +F21*in1+ F22*in2 + u2;
      in = { out0, out1, out2};
    }

    );
  pb.vertices.resize(m.nodes.size() );
  classification classif = {3, 10};
  for( size_t i = 0; i < pb.vertices.size(); ++i){
    pb.vertices[i] = {&classif, i};
  }
  size_t ntets = 0;
  for (int i=0; i < m.elems.size(); ++i){
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
  std::cout << "NTETS " << ntets;
  std::cout << pb.tets[0] << std::endl;

  /*
  pb.coordinates[0][0] = {0.,0.,0.};
  pb.coordinates[0][1] = {1./3.,0.,0.};
  pb.coordinates[0][2] = {0.,1.,0.};
  pb.coordinates[0][3] = {0.,0.,1.};
  
  pb.coordinates[1].resize(4);
  pb.coordinates[1][0] = {0.,0.,0.};
  pb.coordinates[1][1] = {3.,0.,0.};
  pb.coordinates[1][2] = {0.,2.,0.};
  pb.coordinates[1][3] = {0.,0.,2.};
  
  
  pb.vertices ={ { &classif , 1}, { &classif , 2}, { &classif , 3 },  { &classif , 0 }  };
  pb.tets = { {&classif, { &pb.vertices[0],  &pb.vertices[1],  &pb.vertices[2], &pb.vertices[3] } }  };
  */

  
 
 
 

  const Hyperelastic< scal_t > & law = lawneo;
  pb.law = &law;
  
  // int ne = 1 << 26;
  int ne = ntets;

  // pb.tets = std::vector< fem_t::tetrahedre>(ne,  fem_t::tetrahedre{&classif, { &pb.vertices[0],  &pb.vertices[1],  &pb.vertices[2], &pb.vertices[3] } } );
  
  __gnu_parallel::_Parallelism par_tag = __gnu_parallel::_Parallelism::parallel_unbalanced; //130

  pb.setFops();

  // Forcing instanciation of templates
  /*
  fem_t::tensor_t F0={1.,0.,0.,0.,1.,0.,0.,0.,1.};
  auto PKI0= lawneo.PKI(F0);
  auto PKI1= lawst.PKI(F0);
  auto F1 = pb.F(pb.tets[0], 0, 1);
  std::cout << "F elem1 v1 " << std::endl;
  std::cout << F1 << std::endl;
  auto F2 = pb.Fv2(pb.tets[0], 0, 1);
  std::cout << "F elem1 v2 " << std::endl;
  std::cout << F2 << std::endl;
  auto F3 = pb.Fv3(pb.tets[0],  1);
  std::cout << "F elem1 v3 " << std::endl;
  std::cout << F3 << std::endl;
  */
  
  omp_set_num_threads(8);
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
	 //std::cout << Fint[0] << std::endl;
	 //std::cout << Fint[ne-1] << std::endl;
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
  std::cout <<  " parallel transform (base law ) "<< omp_get_max_threads() << " " << time_parallel.count() << std::endl;
  std::cout <<  " serial   transform (base law ) "<< 1  << " " << time_serial.count() << std::endl;
  std::cout <<  " serial   speed up "<<  " " << time_serial.count()/time_parallel.count() << std::endl;
   
 

   
}
