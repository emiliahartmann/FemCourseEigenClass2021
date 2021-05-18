//
//  TopologyTetrahedron.cpp
//  FemSC
//
//  Created by Karolinne Oliveira Coelho on 4/28/18.
//
//

#include <stdio.h>
#include "TopologyTetrahedron.h"
#include "tpanic.h"

// Number of sides associated with tetrahedron elements elements
const int TopologyTetrahedron::nSides;

// Number of corner nodes associated with tetrahedron elements
const int TopologyTetrahedron::nCorners;

// Dimension of tetrahedron elements
const int TopologyTetrahedron::Dimension;


int TopologyTetrahedron::NSideNodes(int side)
{
    if (side>14) {
        std::cout << "TopologyTetrahedron::NSideNodes: Bad parameter side" << std::endl;
        DebugStop();
        return EXIT_FAILURE;
    }
    
    int nsidenodes[15] = {1,1,1,1,2,2,2,2,2,2,3,3,3,3,4};
    return nsidenodes[side];
}

// local node index of a node associated with a side
int TopologyTetrahedron::SideNodeLocIndex(int side, int node)
{
    int FaceNodes[4][3] = {{0,1,2},{0,1,3},{1,2,3},{0,2,3}};
    int SideNodes[6][2] = {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
    if(side<4 && node==0)
        return side;
    else if(side>=4 && side<10 && node<2)
        return SideNodes[side-4][node];
    else if(side>=10 && side<14 && node<3)
        return FaceNodes[side-10][node];
    else if(side==14 && node<4)
        return node;
    
    std::cout << "TopologyTetrahedron::SideNodeIndex inconsistent side or node" << std::endl;
    return EXIT_FAILURE;
}

// return the enumerated element type
MElementType TopologyTetrahedron::Type(){
    return ETetraedro;
}
