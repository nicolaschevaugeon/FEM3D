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
  };

  mesh import_GMSH( std::istream & infile);
  
}//end namespace gmsh_import

#endif
