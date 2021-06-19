//
//  ShapeQuad.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "Shape1d.h"
#include "ShapeQuad.h"

using namespace std;

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeQuad::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){
    
// Implementação da Emilia
    if(xi.size() != Dimension) DebugStop();
    if(orders.size() != nSides) DebugStop();

    int nshape = NShapeFunctions(orders);
    int nsides = nSides;

    orders.resize(nSides);
    phi.resize(nshape);
    dphi.resize(2, nshape);

    VecDouble csi(1);
    csi[0] = xi[0];
    VecDouble eta(1);
    eta[0] = xi[1];

    phi[0] = (1. - csi[0]) * (1. - eta[0])/4.;
    phi[1] = (1. + csi[0]) * (1. - eta[0])/4.;
    phi[2] = (1. + csi[0]) * (1. + eta[0])/4.;
    phi[3] = (1. - csi[0]) * (1. + eta[0])/4.;

    dphi(0,0) = - 1. * (1. - eta[0])/4.;
    dphi(0,1) = + 1. * (1. - eta[0])/4.;
    dphi(0,2) = + 1. * (1. + eta[0])/4.;
    dphi(0,3) = - 1. * (1. + eta[0])/4.;

    dphi(1,0) = - 1. * (1. - csi[0])/4.;
    dphi(1,1) = - 1. * (1. + csi[0])/4.;
    dphi(1,2) = + 1. * (1. + csi[0])/4.;
    dphi(1,3) = + 1. * (1. - csi[0])/4.;    


//  Implementando as funcoes de forma para os lados (manualmente)
//    int ip;
//    for (ip = 4; ip < 8; ip++){
//        if(orders[ip] == 2)
//        {
//        // implementar aqui
//        phi[4] = (1. - csi[0] * csi[0]) * (1. - eta[0])/2.;
//        phi[5] = (1. - eta[0] * eta[0]) * (1. + csi[0])/2.;
//        phi[6] = (1. - csi[0] * csi[0]) * (1. + eta[0])/2.;
//        phi[7] = (1. - eta[0] * eta[0]) * (1. - csi[0])/2.;

//        dphi(0,4) = - csi[0] * (1. - eta[0]);
//        dphi(0,5) = + 1. * (1. - eta[0])/2.;
//        dphi(0,6) = - csi[0] * (1. + eta[0]);
//        dphi(0,7) = - 1. * (1. - eta[0])/2.;

//        dphi(1,4) = - 1. * (1. - csi[0] * csi[0])/2.;
//        dphi(1,5) = - eta[0] * (1. + csi[0]);
//        dphi(1,6) = + 1. * (1. - csi[0] * csi[0])/2.;
//        dphi(1,7) = - eta[0] * (1. - csi[0]);  
//        }
//  
//        else if(orders[ip] != 1) DebugStop();
//        }

//    if(orders[8] == 2) // Caso separado o lado 8
//        {
//        phi[8] = (1. - csi[0] * csi[0]) * (1. - eta[0] * eta[0]);
//        dphi(0,8) = -2. *  csi[0] * (1. - eta[0] * eta[0]);
//        dphi(1,8) = - 2. * eta[0] * (1. - csi[0] * csi[0]);

//        }
//    else if(orders[8] != 1) DebugStop();

// Implementando as funcoes de forma para os lados (automatica)

    int count = 4;
    int ip;
    for (ip = 4; ip < 8; ip++){
        if(orders[ip] == 2)
        {
        int id= ip % 4;
        int id2 = (ip + 1) % 4;
        int id3 = (ip + 2) % 4;
        phi[count] = 4. * phi[id] * (phi[id2] + phi[id3]);
        dphi(0, count) = 4. * (dphi(0, id) * (phi[id2] + phi[id3]) + (dphi(0, id2) + dphi(0, id3)));
        dphi(1, count) = 4. * (dphi(1, id) * (phi[id2] + phi[id3]) + (dphi(1, id2) + dphi(1, id3)));
        count++;
//        std :: cout << "contador" << count << endl;
        }  
        else if(orders[ip] != 1) DebugStop();
    }

    // Caso separado o lado 8
    if(orders[8] == 2)
    {
        phi[count] = 16. * phi[0] * phi[2];
        dphi(0, count) = 16. * (dphi(0, 0) * phi[2] + phi[0] * dphi(0, 2));
        dphi(1, count) = 16. * (dphi(1, 0) * phi[2] + phi[0] * dphi(1, 2));
        count++;
    }
    else if(orders[8] != 1) DebugStop();
    if(count != nshape) DebugStop();

// -----------------------

    for (int i = 0; i < orders.size(); i++)
    {
        if (orders[i] < 0) {
            std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
            DebugStop();
        }
    }
    if (orders[0] > 1 || orders[1] > 1 || orders[2] > 1 || orders[3] > 1) {
        std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }

    auto nf = NShapeFunctions(orders);

    if (orders[nf-1] > 2) {
        std::cout << "ShapeQuad::Shape, only implemented until order = 2" << std::endl;
        DebugStop();
    }

//    std::cout << "Please implement me\n";
//    DebugStop();
}

/// returns the number of shape functions associated with a side
int ShapeQuad::NShapeFunctions(int side, int order){
    if(order < 1 || order >2) DebugStop();
    if(side<4)
        return 1;//0 a 4
    else if(side<8)
        return (order-1);//6 a 14
    else if(side==8)
        return ((order-1)*(order-1));
    
    std::cout << "ShapeQuad::NShapeFunctions, bad parameter side " << side << std::endl;
    DebugStop();
    
    return 0;
}

/// returns the total number of shape functions
int ShapeQuad::NShapeFunctions(VecInt &orders){
    
    int res=4;
    for(int in=4; in<orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }
    
    return res;
}
