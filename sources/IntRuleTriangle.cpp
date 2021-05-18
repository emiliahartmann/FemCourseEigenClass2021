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
            fPoints.resize(1,Dimension());
            fWeights.resize(1);
            fPoints(0,0) = 1./3.;    
            fPoints(0,1) = 1./3;    
            fWeights(0) = 1./2;
            break;
        case 2:
            fPoints.resize(3,Dimension());
            fWeights.resize(3);
//          Primeira opcao (pontos nas arestas)
//            fPoints(0,0) = 1./2;    
//            fPoints(1,0) = 0.;
//            fPoints(2,0) = 1./2;

//            fPoints(0,1) = 1./2;    
//            fPoints(1,1) = 1./2;
//            fPoints(2,1) = 0.;

//            fWeights(0) = 1./6;
//            fWeights(1) = 1./6;
//            fWeights(2) = 1./6;

//          Segunda opcao (pontos dentro do elemento, nas bissetrizes)
            fPoints(0,0) = 1./6;    
            fPoints(1,0) = 2./3;
            fPoints(2,0) = 1./6;

            fPoints(0,1) = 1./6;    
            fPoints(1,1) = 1./6;
            fPoints(2,1) = 2./3;

            fWeights(0) = 1./6;
            fWeights(1) = 1./6;
            fWeights(2) = 1./6;
            break;
        case 3:
            fPoints.resize(4,Dimension());
            fWeights.resize(4);
            fPoints(0,0) = 1./3;    
            fPoints(1,0) = 1./5;
            fPoints(2,0) = 1./5;
            fPoints(3,0) = 3./5;

            fPoints(0,1) = 1./3;    
            fPoints(1,1) = 1./5;
            fPoints(2,1) = 3./5;
            fPoints(3,1) = 1./5;

            fWeights(0) = -27./96;
            fWeights(1) = 25./96;
            fWeights(2) = 25./96;
            fWeights(3) = 25./96;
            break;
        case 4:  
        case 5:
            fPoints.resize(6,Dimension());
            fWeights.resize(6);
            fPoints(0,0) = 949./1440;    
            fPoints(1,0) = 949./1440;
            fPoints(2,0) = 613./2643;
            fPoints(3,0) = 613./2643;
            fPoints(4,0) = 573./5255;
            fPoints(5,0) = 573./5255;

            fPoints(0,1) = 613./2643;    
            fPoints(1,1) = 573./5255;
            fPoints(2,1) = 949./1440;
            fPoints(3,1) = 573./5255;
            fPoints(4,0) = 949./1440;
            fPoints(5,0) = 613./2643;

            fWeights(0) = 1./6;
            fWeights(1) = 1./6;
            fWeights(2) = 1./6;
            fWeights(3) = 1./6;
            fWeights(4) = 1./6;
            fWeights(5) = 1./6;
            break;

        break;
        default:
            DebugStop();
            break;
        }

}
