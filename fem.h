#ifndef _fem_
#define _fem_

#include <vector>
#include "tensor2_9cm.h"
#include "hyperelastic.h"

class classification{
public :
  size_t dim;
  size_t id;
};
using h_classification = classification *;




template <class Scal >
using h_coordinate = coordinate<Scal> *;


class meshentity{
public :
  h_classification h_class;
};

class vertex : public meshentity {
public:
  vertex()  =default;
 vertex( h_classification classif, size_t _index, size_t _base_mesh_id ): meshentity{classif}, index{_index}, base_mesh_id(_base_mesh_id){} 
  size_t index;
  size_t base_mesh_id;
};



using h_vertex = vertex *;

template < class Scal >
class basis : public std::array< Scal , 9>{
public: 
  basis( const coordinate<Scal > &X0, const coordinate<Scal > &X1, const coordinate<Scal > &X2):
  std::array<Scal, 9>{ {X0[0], X0[1], X0[2], X1[0], X1[1], X1[2], X2[0], X2[1], X2[2]  } } {}
};

/*

field{
  space;
  map <entity, coefficients>;
}

vector3d_field_coefficient{
std::array<int, 3>     indexes_in_equation_array;
std::array< size_t, 3> indexes_in_unknown array;
unordered_map <dof_id, std::pair< entity *, dir > >
}
*/

/*
class vector3d_field_coefficient{
  enum type {free, fixed}
  size_t number;
  std::unordered_map
 
  }*/

//tensor2_9cm<double > 

template <class Scal >
class fem{
public:
  using scal_t =Scal;
  using coord_t = coordinate<Scal>;
  using basis_t = basis<Scal>;
  using tensor_t =  tensor2_9cm<Scal >;
  using mat34cm_t = mat34cm<Scal>;
  class tetrahedre : public meshentity {
   
  public:
    tetrahedre (h_classification classif, const  std::array< h_vertex, 4> & _vertices, size_t _id):  meshentity{classif}, // the following does not work on g++ 4.8.2 ... vertices{ _vertices },
id(_id) {
  //std::copy(_vertices.begin(), _vertices.begin(), vertices.begin() );
  vertices[0] = _vertices[0];
  vertices[1] = _vertices[1];
  vertices[2] = _vertices[2];
  vertices[3] = _vertices[3];

    }
    std::array< h_vertex, 4> vertices;
    size_t id;
  };

  coord_t get_coordinate(const vertex &v, size_t  confindex){
    return coordinates[confindex][v.index];
  };
  
  std::array< coord_t, 4> get_coordinates(const tetrahedre& tet,  size_t confindex){
    return std::array<coord_t, 4> {get_coordinate(*(tet.vertices[0]), confindex), get_coordinate( *(tet.vertices[1]), confindex), get_coordinate(*(tet.vertices[2]), confindex), get_coordinate(*(tet.vertices[3]), confindex) };
  }

  mat34cm_t get_coordinates34(const tetrahedre & tet, size_t confindex){
    const coord_t & x0 = get_coordinate(*(tet.vertices[0]), confindex);
    const coord_t & x1 = get_coordinate(*(tet.vertices[1]), confindex);
    const coord_t & x2 = get_coordinate(*(tet.vertices[2]), confindex);
    const coord_t & x3 = get_coordinate(*(tet.vertices[3]), confindex);
    return mat34cm_t
      {   x0[0], x0[1], x0[2],
	  x1[0], x1[1], x1[2],
	  x2[0], x2[1], x2[2],
	  x3[0], x3[1], x3[2]};
  }
  
  //std::vector<basis  >       get_basis( size_t confindex, const   std::vector< coordinate > & );
  basis_t get_basis(const tetrahedre &tet, size_t confindex ){
    std::array< coord_t, 4> coords = get_coordinates(tet,  confindex);
    return basis_t { coords[1]-coords[0], coords[2] -coords[0], coords[3] - coords[0]  };
  }
  
  //Calcul F entre la conf 0 et 1
  tensor_t F (const tetrahedre &tet,  size_t confindex0, size_t confindex1 ){
    
    tensor_t basis1 = { get_basis(tet, confindex1) };
    tensor_t basis0 = { get_basis(tet, confindex0) };
    return basis1*invert(basis0);
  }

  void setFops(){
    Fops.resize(tets.size());
    std::transform(tets.begin(), tets.end(), Fops.begin(), [this]( const tetrahedre &tet){  return this-> getFop(tet, 0);} );
  }
  
  mat34cm_t getFop(const tetrahedre &tet, size_t confindex0  ){
    tensor_t dualbasis  = invert_transpose (tensor_t { get_basis(tet, confindex0) });
    //tensor_t dualbasis  = invert (tensor_t { get_basis(tet, confindex0) });
    // dualbasis = {1.,0.,0.,0.,1.,0.,0.,0.,1};
    const auto &dbasis = dualbasis.data;
    //std::cout << "dual basis" << std::endl;
    //std::cout << dualbasis << std::endl;
    return mat34cm_t {
      -dbasis[0] - dbasis[3] - dbasis[6],
      -dbasis[1] - dbasis[4] - dbasis[7],
      -dbasis[2] - dbasis[5] - dbasis[8],
	dbasis[0],
	dbasis[1],
	dbasis[2],
	dbasis[3],
	dbasis[4],
	dbasis[5],
	dbasis[6],
	dbasis[7],
	dbasis[8]
	};
    // std::cout << "Fop" << std::endl;
    //std::cout << res << std::
    //    return res;
  }

  tensor_t Fv2(const tetrahedre &tet, size_t confindex0,  size_t confindex1  ){
    const auto x = get_coordinates34(tet, confindex1);
    const auto Fop = getFop(tet,  confindex0  );
    return apply( Fop, x);
  }

  tensor_t Fv3(const tetrahedre &tet, size_t confindex1  ){
    const auto x = get_coordinates34(tet, confindex1);
    return apply( Fops[tet.id], x);
  }

  mat34cm<scal_t > internal_forces(const tetrahedre &tet, size_t confindex1  ){
    const auto xt = get_coordinates34(tet, confindex1);
    const auto &Fopst = Fops[tet.id];
    const auto Ft = apply( Fopst, xt);
    const auto Pt = law->PKI(Ft);
    return apply( Fopst, Pt);
  }
  //std::vector < vertex >
  std::array < std::vector < coord_t > , 2 > coordinates;
  std::vector< vertex > vertices;
  std::vector < tetrahedre > tets;
  std::vector < mat34cm_t >  Fops;
  const Hyperelastic< scal_t> *law = nullptr;
  //std::vector < edge >
  //std::vector < triangle >
  
};

//template <class Scal> 
std::ostream & operator << (std::ostream & out, const typename fem<double>::tetrahedre &t){
  out << t.vertices[0]<< " "  << t.vertices[1] << " " <<t.vertices[2] << " "<< t.vertices[3];
  return out; 
}

std::ostream & operator << (std::ostream & out, const typename fem<float>::tetrahedre &t){
  out << t.vertices[0]<< " "  << t.vertices[1] << " " <<t.vertices[2] << " "<< t.vertices[3];
  return out; 
}




#endif
