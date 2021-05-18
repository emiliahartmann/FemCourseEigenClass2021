//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#ifndef Geom1d_h
#define Geom1d_h

#include "Topology1d.h"
#include "GeoElementSide.h"

/**
 @brief One dimensional geometry
 @ingroup geometrymap
 */
class Geom1d : public Topology1d
{
public:
    
    // Constructor
    Geom1d();
    
    // Destructor
    ~Geom1d();
    
    // Copy constructor
    Geom1d(const Geom1d &copy);
    
    // Operator of copy
    Geom1d &operator=(const Geom1d &copy);
    
    // Computes the shape functions associated with the geometric map
    static void Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi);
    
    // Computes the value of x for a given point in parameter space as a function of corner coordinates
    static void X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x);
    
    // Computes the value of x and gradx for a given point in parameter space
    static void GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx);
    
    // Return the number of nodes of the template
    static int NumNodes();
    
    // Set the node indices of the element
    void SetNodes(const VecInt &nodes);
    
    // Set the node indices of the element
    void GetNodes(VecInt &nodes) const;
    
    // Return the index of a node
    int NodeIndex(int node) const;
    
    // Return the neighbour along side
    GeoElementSide Neighbour(int side) const;
    
    // Initialize the neighbour data structure
    void SetNeighbour(int side, const GeoElementSide &neighbour);
    
protected:
    
    // Indexes of the nodes of the geometry
	VecInt fNodeIndices;
    
    // Vector of neighbours
    GeoElementSide fNeighbours[nSides];
};

#endif /* Geom1d_h */
