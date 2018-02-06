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
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

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
  using scal_t = double;
  if (argc < 2){
    std::cerr  << "Usage : " << argv[0] << " meshfilename.msh [numthreads]"<< std::endl;
    return 1;
  }

  
  int numthreads = 4 ;
  if (argc == 3) numthreads = std::stoi (argv[2]);
  omp_set_num_threads(numthreads);
 
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
    pb.vertices[i] = {&classif, i, i+1};
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
  
  // setting deformed configuration  nodes coordinates
  for_each( pb.coordinates[1].begin(), pb.coordinates[1].end(), []( coordinate< scal_t> &X  ){
      const coordinate<scal_t> u ={2.,3.,4.};
      const tensor2_9cm< scal_t> F = {1.3, 0.5, 0.2,1.2, 1.01, 0.6, 2., 0.7, 0.37 };
      X =   F*X +u;
     }
    
    );
  /*
  fs::copy_file(meshfilename, "result.msh", fs::copy_options::overwrite_existing );
  append_vector_view_gmsh( "result.msh", "displacement", 0.,0, [&pb] (const vertex &v) {  const auto  u = pb.coordinates[1][v.index] -pb.coordinates[0][v.index];         return std::array<float, 3>{ (float)u.data[0], (float)u.data[1], (float)u.data[2]};
    },
    [&pb] (const vertex &v) { return v.base_mesh_id;},
    pb.vertices.begin(), pb.vertices.end()
    );
  */
   
  std::cout << "NTETS " << ntets;
  std::cout << pb.tets[0] << std::endl;
 
  //Hyperelastic_StVenantKirchhoff < scal_t > lawst = {1.d,2.d};
  Hyperelastic_NeoHookeen < scal_t > lawneo = {1.d,2.d};
  
 

  const Hyperelastic< scal_t > & law = lawneo;
  pb.law = &law;
  
  int ne = ntets;
 
  __gnu_parallel::_Parallelism par_tag = __gnu_parallel::_Parallelism::parallel_unbalanced; //130

  pb.setFops();
 
  auto now = std::chrono::high_resolution_clock::now;
  typedef std::chrono::duration<double, std::milli> dt_t;

  {
    std::vector< mat34cm< scal_t > >  Fint(ne);
    auto t_start = now();
    __gnu_parallel::transform(pb.tets.begin(), pb.tets.end(), Fint.begin(),  [&pb](const fem_t::tetrahedre & tet ){
	return pb.internal_forces( tet, 1);
      },
      par_tag);
    auto t_end = now();
    dt_t time_parallel = t_end-t_start;
    std::cout <<  " parallel transform (base law ) "<< omp_get_max_threads() << " " << time_parallel.count() << std::endl;
    std::cout << Fint[0] << std::endl;
    std::cout << Fint[ne-1] << std::endl;
  }

  {
    std::vector < scal_t>  Fint(pb.vertices.size()*3, 0. );
    auto t_start = now();
    std::for_each(pb.tets.begin(), pb.tets.end(),  [&pb, &Fint](const fem_t::tetrahedre & tet ){
	const mat34cm< scal_t >	Fint_e  = pb.internal_forces( tet, 1);
	for (int i=0; i < 4; ++i){
	  const size_t index = tet.vertices[i]->index;
	  Fint[3* index]   +=  Fint_e.data[i*3 ];
	  Fint[3* index+1] +=  Fint_e.data[i*3+1 ];
	  Fint[3* index+2] +=  Fint_e.data[i*3+2 ];
	}
      });
    auto t_end = now();
    dt_t time_parallel = t_end-t_start;
    std::cout <<  " assemble Fint serial " << " " << time_parallel.count() << std::endl;
  }

   {
     using F_type = std::vector < scal_t>;
     F_type Fint(pb.vertices.size()*3, 0. );
     const size_t nth = omp_get_max_threads();
     std::vector< F_type>  Fint_th (  nth, F_type(pb.vertices.size()*3, 0.) );
     auto t_start = now();
     std::vector < ITER1> b_t(nth);
     std::vector < ITER1> e_t(nth);
     const size_t l = std::distance(b1,e1);
     for (size_t i = 0; i < nth; ++i){
       b_t[i] = b+i*l/nth;
       e_t[i] = b+(i+1)*l/nth;
     }
     e_t[nth-1] = e;
     
#omp parallel
     {
    
#pragma omp parallel for
  for (size_t i =0; i< nth; ++i)
    std::transform(b1_t[i], e1_t[i],  b2_t[i], op  );
       
       std::for_each(pb.tets.begin(), pb.tets.end(),  [&pb, &Fint](const fem_t::tetrahedre & tet ){
	   const mat34cm< scal_t >	Fint_e  = pb.internal_forces( tet, 1);
	   for (int i=0; i < 4; ++i){
	     const size_t index = tet.vertices[i]->index;
	     Fint[3* index]   +=  Fint_e.data[i*3 ];
	     Fint[3* index+1] +=  Fint_e.data[i*3+1 ];
	     Fint[3* index+2] +=  Fint_e.data[i*3+2 ];
	   }
	 });
       
       dt_t time_parallel = t_end-t_start;
       std::cout <<  " assemble Fint serial " << " " << time_parallel.count() << std::endl;
     }
     
   }
	 
 


  
 
  return 0;
   
}
