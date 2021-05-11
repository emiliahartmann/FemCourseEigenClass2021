/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream> 
#include "IntRuleTriangle.h"
#include "tpanic.h"

IntRuleTriangle::IntRuleTriangle(){

}

IntRuleTriangle::IntRuleTriangle(int order) {
    SetOrder(order);
//    DebugStop();
}

void IntRuleTriangle::SetOrder(int order) {
    fOrder=order;
    if(order < 0 || order > MaxOrder());
    switch (order){
        case 0:
        case 1:
//            fPoints.resize(1,Dimension());
//            fWeights.resize(1);
//            fPoints(0,0) =  1./3; 
//            fPoints(0,1) =  1./3;
//            fWeights(0) = 1.;
            fPoints.resize(3,Dimension());
            fWeights.resize(3);
            fPoints(0,0) =  -1./3.;    
            fPoints(1,0) =  1./3;
            fPoints(0,1) =  0.;    
            fPoints(1,1) =  0.;
            fWeights(0) = 1./3;
            fWeights(1) = 1./3;
            DebugStop();
            break;
        case 2:
        case 3:
            fPoints.resize(3,Dimension());
            fWeights.resize(3);
            fPoints(0,0) =  0.;    
            fPoints(1,0) =  0.;
            fPoints(2,0) =  1.;
            fPoints(0,1) =  0.;    
            fPoints(1,1) =  0.;
            fPoints(2,1) =  1.;
            fWeights(0) = 1./3;
            fWeights(1) = 1./3;
            fWeights(2) = 1./3;
            break;
        default:
            DebugStop();
            break;
        }

    DebugStop();
}
