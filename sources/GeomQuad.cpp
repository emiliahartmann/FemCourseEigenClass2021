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
    if(xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) DebugStop();

    double csi = xi[0];
    double eta = xi[1];

    phi[0] = (1. - csi) * (1. - eta)/4.;
    phi[1] = (1. + csi) * (1. - eta)/4.;
    phi[2] = (1. + csi) * (1. + eta)/4.;
    phi[3] = (1. - csi) * (1. + eta)/4.;

    dphi(0,0) = - 1. * (1. - eta)/4.;
    dphi(0,1) = + 1. * (1. - eta)/4.;
    dphi(0,2) = + 1. * (1. - eta)/4.;
    dphi(0,3) = - 1. * (1. - eta)/4.;
    dphi(1,0) = - 1. * (1. - csi)/4.;
    dphi(1,1) = - 1. * (1. + csi)/4.;
    dphi(1,2) = + 1. * (1. + csi)/4.;
    dphi(1,3) = + 1. * (1. - csi)/4.;

//    DebugStop();
}

void GeomQuad::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() != NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();

    int nrow=NodeCo.rows();
    int ncol=NodeCo.cols();

    x.resize(nrow);
    x.setZero();

    VecDouble phi(4);
    MatrixDouble dphi(Dimension,4);
    Shape(xi, phi, dphi);
    for(int i = 0; i < ncol; i++){
        for(int j = 0; j < nrow; j++){
            x[j] += NodeCo(j,i) * phi[i];      
        } 
    }

//    DebugStop();
}

void GeomQuad::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {

    if(xi.size() != Dimension) DebugStop();
    if(x.size() != NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();

    int nrow=NodeCo.rows();
    int ncol=NodeCo.cols();

    gradx.resize(nrow,2);
    gradx.setZero();
    x.resize(nrow);
    x.setZero();

    VecDouble phi(4);
    MatrixDouble dphi(Dimension,4);
    Shape(xi, phi, dphi);
    for(int i = 0; i < ncol; i++){
        for(int j = 0; j < nrow; j++){
            x[j] += NodeCo(j,i) * phi[i];
            gradx(j, 0) += NodeCo(j, i) * dphi(0, i);       
        } 
    }

//    DebugStop();
}

void GeomQuad::SetNodes(const VecInt &nodes) {
    if(nodes.size() != nCorners) DebugStop();
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
