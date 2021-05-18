//
//  ReadGmsh.h
//  FemSC
//
//  Created by Philippe Devloo on 17/04/18.
//

#ifndef ReadGmsh_h
#define ReadGmsh_h

#include <string>
#include <vector>
#include <map>
#include "GeoMesh.h"

/**
@brief Read the definition of a geometry from a Gmsh file
@ingroup geometry
*/
class ReadGmsh
{
    
public:

    // Constructor of ReadGmsh
    ReadGmsh() : fMaterialDataVec(3)
    {
        
    }
    // Reads the mesh contained in the file and fill the geometric mesh
    void Read(GeoMesh &gmesh, const std::string &filename);
    
protected:
    
    /** @brief MaterialVec */
    /** Structure of both, physical entities dimension and names */
    std::vector<std::map<int,std::string> > fMaterialDataVec;

};
#endif /* ReadGmsh_h */
