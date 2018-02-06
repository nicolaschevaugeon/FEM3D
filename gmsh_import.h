#ifndef _gmsh_import_
#define _gmsh_import_
#include <vector>
#include <iostream>

namespace gmsh_import{
  struct node{
    size_t gmshid;
    double x, y, z;
  };

  struct element{
    int number;
    int type;
    std::vector<int > tags; //geo_class is tags[1]
    std::vector<int > gmshnodeids;
  };

  struct mesh{
    std::vector< node > nodes;
    std::vector< element > elems;
    void clear(){
      nodes.clear();
      elems.clear();
    }
  };

  mesh import_GMSH( std::istream & infile);
  
}//end namespace gmsh_import


#include <functional>
#include <fstream>
#include <algorithm>

template <class VERTEXITERATOR  >
void append_vector_view_gmsh(const std::string & gmshfilename, const std::string & viewname, float timevalue,  int timestep, std::function< std::array<float, 3 >   ( const  typename VERTEXITERATOR::value_type  & v ) >   vvalue, std::function <size_t  ( const  typename VERTEXITERATOR::value_type  & v ) >   vid, const VERTEXITERATOR begin, const VERTEXITERATOR end   )
{
  std::ifstream in(gmshfilename.c_str() );
  if(!in.is_open()) {
    std::cout << "file " << gmshfilename << " Does not exist " << std::endl;
    return;
  }
  in.close();
  std::ofstream out( gmshfilename.c_str(), std::ofstream::app);
  using vertex_type= typename VERTEXITERATOR::value_type ;
  out << "$NodeData" << std::endl;
  out << "1" << std::endl;
  out << "\""<< viewname << "\"" << std::endl;
  out <<"1" << std::endl;
  out <<timevalue << std::endl;
  out<< "3" << std::endl; //nb int tag
  out << timestep << std::endl; //time step
  out << "3"<< std::endl;  //scal (3 for vector)
  out << std::distance(begin, end) << std::endl; //nb values
  std::for_each (begin, end, [&vvalue, &vid, &out](vertex_type & v) {
    auto val = vvalue(v);
    out << vid(v) << " " << val[0] << " " << val[1]<< " " << val[2] << std::endl; });
  out << "$EndNodeData" << std::endl;
}

#endif
