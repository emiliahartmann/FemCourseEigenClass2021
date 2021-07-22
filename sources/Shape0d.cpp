/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Shape0d.h"
//#include "TMatrix.h"

void Shape0d::Shape(const VecDouble& xi, VecInt& orders, VecDouble& phi, MatrixDouble& dphi) {
    // if(xi.size() != Dimension) DebugStop();
    // if(orders.size() != nSides) DebugStop();

    int nshape = NShapeFunctions(orders);
    // int nsides = nSides;

    // phi.resize(nshape);
    // dphi.resize(1,nshape);
    // orders.resize(nSides); 

    phi[0] = 1.;
    // dphi(0,0) = 0.;
}

int Shape0d::NShapeFunctions(int side, int order) {

    if (side == 0) {
        // if(order != 1) DebugStop(); // descomentar isso depois
        return 1;
    }
    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return -1;
}

int Shape0d::NShapeFunctions(VecInt &orders) {
    int n = orders.size();
    if(n != 1) DebugStop();
    int nshape = 0;
    int i = 0;
    for (i = 0; i < n; i++) {
        nshape = nshape + NShapeFunctions(i, orders[i]);
    }
    return nshape;
}
