/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "GeomQuad.h"

GeomQuad::GeomQuad() {
}

GeomQuad::~GeomQuad() {
}

GeomQuad::GeomQuad(const GeomQuad &copy) {
    fNodeIndices = copy.fNodeIndices;
}

GeomQuad& GeomQuad::operator=(const GeomQuad& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void GeomQuad::Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi) {
    // if(xi.size() <=0 || xi.size() > Dimension) DebugStop();
   if(xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) DebugStop();

    // phi.resize(nCorners);
    // dphi.resize(Dimension, nCorners);

    double csi = xi[0];
    double eta = xi[1];

    phi[0] = (1. - csi) * (1. - eta)/4.;
    phi[1] = (1. + csi) * (1. - eta)/4.;
    phi[2] = (1. + csi) * (1. + eta)/4.;
    phi[3] = (1. - csi) * (1. + eta)/4.;

    dphi(0,0) = - 1. * (1. - eta)/4.;
    dphi(0,1) = + 1. * (1. - eta)/4.;
    dphi(0,2) = + 1. * (1. + eta)/4.;
    dphi(0,3) = - 1. * (1. + eta)/4.;

    dphi(1,0) = - 1. * (1. - csi)/4.;
    dphi(1,1) = - 1. * (1. + csi)/4.;
    dphi(1,2) = + 1. * (1. + csi)/4.;
    dphi(1,3) = + 1. * (1. - csi)/4.;

//    DebugStop();
}

void GeomQuad::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() < NodeCo.rows()) DebugStop(); // x.size() != NodeCo.rows()
    if(NodeCo.cols() != nCorners) DebugStop();

    // int nrow=NodeCo.rows();
    // int ncol=NodeCo.cols();

    // if (x.size() < nrow){
    //     x.resize(nrow);
    // }
    
    // x.setZero();

    VecDouble phi(nCorners);  // nCorners == 4
    MatrixDouble dphi(Dimension, nCorners); // nCorners == 4

    Shape(xi, phi, dphi);
    // for(int i = 0; i < nCorners; i++){
    //     for(int j = 0; j < Dimension; j++){
    //         x[j] += NodeCo(j,i) * phi[i];      
    //     }         
    for(int i = 0; i < nCorners; i++){
        for(int j = 0; j < Dimension; j++){
            x[j] += NodeCo(j,i) * phi[i];      
        }     
    // for(int i = 0; i < ncol; i++){
    //     for(int j = 0; j < nrow; j++){
    //         x[j] += NodeCo(j,i) * phi[i];      
    //     } 
    }

//    DebugStop();
}

void GeomQuad::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() < NodeCo.rows()) DebugStop();  
    if(NodeCo.cols() != nCorners) DebugStop();

    int nrow=NodeCo.rows();
    int ncol=NodeCo.cols();
    
    gradx.resize(nrow, Dimension);
    // if (gradx.cols()<nrow) gradx.resize(nrow,1);
    // if (gradx.cols() < nrow) gradx.resize(nrow,Dimension);
    x.resize(nrow);
    gradx.setZero();
    x.setZero();


   
    VecDouble phi(nCorners);
    MatrixDouble dphi(Dimension, nCorners);
    Shape(xi, phi, dphi);

    for(int i = 0; i < nCorners; i++){
        for(int j = 0; j < Dimension; j++){
            x[j] += NodeCo(j,i) * phi[i];
            gradx(j, 0) += NodeCo(j, i) * dphi(0, i);
            gradx(j, 1) += NodeCo(j, i) * dphi(1, i);       
        } 
    } // nrow ou Dimension ?
}

void GeomQuad::SetNodes(const VecInt &nodes) {
    if(nodes.size() != nCorners) {DebugStop();}
    fNodeIndices = nodes;
}

void GeomQuad::GetNodes(VecInt &nodes) const{
    nodes = fNodeIndices;
}

int GeomQuad::NodeIndex(int node) const {
    return fNodeIndices[node];
}

int GeomQuad::NumNodes() {
    return nCorners;
}

GeoElementSide GeomQuad::Neighbour(int side) const {
    return fNeighbours[side];
}

void GeomQuad::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
