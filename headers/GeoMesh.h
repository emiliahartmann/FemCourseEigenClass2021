//
//  GeoMesh.h
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#ifndef GeoMesh_h
#define GeoMesh_h

#include "GeoNode.h"
#include <string>
class CompMesh;
class GeoElement;

/**
@brief Geometric mesh class
 
 This class holds a collection of geometric elements and nodes
@ingroup geometry
*/
class GeoMesh
{
    
    
    // Vector of nodes
    std::vector<GeoNode> Nodes;
    
    // Vector of element pointers
    std::vector<GeoElement *> Elements;
    
    // Computational mesh associated
    CompMesh *Reference;
    
    // Dimension of the geometric domain
    int fDim;
    
public:
    
    // Default Constructor of GeoMesh
    GeoMesh();
    
    // Constructor of GeoMesh
    GeoMesh(const GeoMesh &);
    
    // Operator of copy 
    GeoMesh &operator=(const GeoMesh &);
    
    // Set the number of nodes on the grid
    void SetNumNodes(int nnodes);
    
    // Set the number of elements on the grid
    void SetNumElements(int numelements);
    
    // Number of nodes of the mesh
    int NumNodes() const;
    
    // Return the number of elements
    int NumElements() const;
    
    // Return the node associated with a index
    const GeoNode &Node(int node) const;
    
    // Return the node associated with a index
    GeoNode &Node(int node);
    
    // Set computational element pointer
    void SetElement(int elindex, GeoElement *gel);
    
    // Set Dimension
    void SetDimension(int dim){fDim = dim;}
    
    // Get Dimension.
    int Dimension() const {return fDim;}

    // Return the elements associated with a index
    GeoElement *Element(int elindex) const;
    
    // Build the connectivity of the grid
    void BuildConnectivity();
    
    // Sets the reference of the geometric grid to ref
    void SetReference(CompMesh *ref)
    {
        Reference = ref;
    }
    
    // Returns the currently loaded computational grid
    CompMesh *GetReference() const {return Reference;}
    
    // Function to print results
    void Print(std::ostream &out) const;
    
};
#endif /* GeoMesh_h */
