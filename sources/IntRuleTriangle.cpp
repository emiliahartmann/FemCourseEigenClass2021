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
    DebugStop();
}

void IntRuleTriangle::SetOrder(int order) {
    fOrder=order;
    if(order < 0 || order > MaxOrder());
    switch (order){
        case 0:
            DebugStop();
            break;
        case 1:
            DebugStop();
            break;
        case 2:
            DebugStop();
            break;
        default:
            DebugStop();
            break;
        }

    DebugStop();
}
