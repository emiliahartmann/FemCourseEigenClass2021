//
//  VTKGeoMesh.h
//  FemSC
//
//  Created by Philippe Devloo on 17/04/18.
//

#ifndef VTKGeoMesh_h
#define VTKGeoMesh_h

#include <string>
#include "DataTypes.h"

class GeoMesh;
class CompMesh;
class PostProcess;

class VTKGeoMesh
{
    public:
    /** @brief Generate an output of all geomesh to VTK */
    static void PrintGMeshVTK(GeoMesh *gmesh, const std::string &filename);

    /// Generate an output file for the solution and its gradient
    static void PrintCMeshVTK(CompMesh *cmesh, int dim, const std::string &filename);
    
    /// Generate an output file for the solution and its gradient
    static void PrintSolVTK(CompMesh *cmesh, PostProcess &defPostProc, const std::string &filename);
    
};

#endif /* VTKGeoMesh_h */