//
//  Topology1d.cpp
//  FemSC
//
//  Created by Karolinne Oliveira Coelho on 4/28/18.
//
//

#include "Topology1d.h"
#include "tpanic.h"

// Number of sides associated with 1D elements
const int Topology1d::nSides;

// Number of corner nodes associated with 1D elements
const int Topology1d::nCorners;

// Dimension of 1D elements
const int Topology1d::Dimension;


int Topology1d::NSideNodes(int side)
{
    if (side>2) {
        std::cout << "TopologyQuad::NSideNodes: Bad parameter side" << std::endl;
        DebugStop();
        return EXIT_FAILURE;
    }
    
    int nsidenodes[3] = {1,1,2};
    return nsidenodes[side];
}

// local node index of a node associated with a side
int Topology1d::SideNodeLocIndex(int side, int node)
{
    if(side <2 && node == 0)
        return side;
    if(side == 2 && node <2)
        return node;
    std::cout << "Topology1d::SideNodeIndex inconsistent side or node" << std::endl;
    return -1;
}

// return the enumerated element type
MElementType Topology1d::Type(){
    return EOned;
}
