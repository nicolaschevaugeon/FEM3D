#include "gmsh_import.h"
#include <string>


namespace gmsh_import{
  mesh import_GMSH( std::istream & infile){
    mesh themesh;
    std::vector< node > &nodes = themesh.nodes;
    std::vector< element > &elems = themesh.elems;
    std::string line;
    
    // reading nodes
    int nnodes;
    do{
      line.clear();
      std::getline(infile, line);
    }while(infile.good() && ((line.find("$Nodes") == std::string::npos) && (line.find("$NOD") == std::string::npos)   ));
    if (!infile.good() ){
      std::cout << "Can't find $Nodes or $NOD in stream infile " << __FILE__ <<":"<< __LINE__ << std::endl;
      throw;
    }
    infile >> nnodes;
    //std::cout << nnodes << std::endl;
    nodes.reserve( nnodes);
    for (int i=0; i < nnodes; ++i){
      node vk;
      infile >> vk.gmshid >> vk.x >> vk.y >> vk.z;
      nodes.push_back(std::move(vk));
    }
    
    //reading elements
    int nelem;
    do{
      line.clear();
      std::getline(infile, line);
    }while(infile.good() && (line.find("$Elements") == std::string::npos)  &&  (line.find("$ELM") == std::string::npos));
    if (!infile.good() ){
      std::cout << "Can't find $Elements or $ELM in stream infile " << __FILE__ <<":"<< __LINE__ << std::endl;
      throw;
    }
    infile>> nelem;
    //std::cout << nelem << std::endl;
    elems.reserve(nelem);
    for (int i=0; i < nelem; ++i){
      element ek;
      int nbtags;
      
      infile >> ek.number >> ek.type >> nbtags;
      ek.tags.resize(nbtags);
      for (int k =0; k < nbtags; ++k)
	infile >> ek.tags[k];
      int nnodese;
      switch(ek.type){
      case 15:{ nnodese =1; break;}
      case 1 :{ nnodese =2; break;} 
      case 2 :{ nnodese =3; break;} 
      case 4 :{ nnodese =4; break;}   
      default :{
	std::cout << "Error in importGMSH " << __FILE__<< ":"<< __LINE__ << " Gmsh element type "<< ek.type << " unknown" << std::endl; 
	throw;
      }
      }
      ek.gmshnodeids.resize(nnodese);
      for (int m =0; m < nnodese; ++m)
	infile >> ek.gmshnodeids[m];
      elems.push_back(std::move(ek));
    }
    return themesh;
  }

}//end namspace gmsh_import
 
