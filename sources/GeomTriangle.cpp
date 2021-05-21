/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeomTriangle.h"

GeomTriangle::GeomTriangle() {
}

GeomTriangle::~GeomTriangle() {
}

GeomTriangle::GeomTriangle(const GeomTriangle &copy) {
    fNodeIndices = copy.fNodeIndices;

}

GeomTriangle& GeomTriangle::operator=(const GeomTriangle& copy) {
    fNodeIndices = copy.fNodeIndices;

    return *this;
}

void GeomTriangle::Shape(const VecDouble& xi, VecDouble& phi, MatrixDouble& dphi) {
    if(xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) DebugStop();

    double csi = xi[0];
    double eta = xi[1];

    phi[0] = 1. - csi - eta;
    phi[1] = csi;
    phi[2] = eta;

    dphi(0,0) = - 1;
    dphi(0,1) = + 1.;
    dphi(0,2) = 0.;

    dphi(1,0) = - 1.;
    dphi(1,1) = 0.;
    dphi(1,2) = + 1. ;

//    DebugStop();
}

void GeomTriangle::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() != NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();

    int nrow=NodeCo.rows();
    int ncol=NodeCo.cols();

    x.resize(nrow);
    x.setZero();

    VecDouble phi(3);
    MatrixDouble dphi(Dimension,3);
    Shape(xi, phi, dphi);
    for(int i = 0; i < ncol; i++){
        for(int j = 0; j < nrow; j++){
            x[j] += NodeCo(j,i) * phi[i];      
        } 
    }

//    DebugStop();
}

void GeomTriangle::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() != NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();

    int nrow=NodeCo.rows();
    int ncol=NodeCo.cols();

    gradx.resize(nrow,2);
    gradx.setZero();
    x.resize(nrow);
    x.setZero();

    VecDouble phi(3);
    MatrixDouble dphi(Dimension,3);
    Shape(xi, phi, dphi);
    for(int i = 0; i < ncol; i++){
        for(int j = 0; j < nrow; j++){
            x[j] += NodeCo(j,i) * phi[i];
            gradx(j, 0) += NodeCo(j, i) * dphi(0, i);       
        } 
    }

//    DebugStop();
}

void GeomTriangle::SetNodes(const VecInt &nodes) {
    if(nodes.size() != nCorners) DebugStop();
    fNodeIndices = nodes;
}

void GeomTriangle::GetNodes(VecInt &nodes) const  {
    nodes = fNodeIndices;
}

int GeomTriangle::NodeIndex(int node) const {
    return fNodeIndices[node];
}

int GeomTriangle::NumNodes() {
    return nCorners;
}

GeoElementSide GeomTriangle::Neighbour(int side)  const {
    return fNeighbours[side];
}

void GeomTriangle::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
