/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Geom1d.h"

Geom1d::Geom1d() {
}

Geom1d::~Geom1d() {
}

Geom1d::Geom1d(const Geom1d &copy) {
    fNodeIndices = copy.fNodeIndices;
}

Geom1d& Geom1d::operator=(const Geom1d& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void Geom1d::Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi) {
//    std::cout << __PRETTY_FUNCTION__ << std::endl;
    if(xi.size() <=0 || xi.size() > Dimension) DebugStop();
//    if(xi.size() !=Dimension ) DebugStop();
//|| phi.size() !=nCorners || dphi.rows() !=Dimension || dphi.cols() !=nCorners
    phi.resize(nCorners);
    dphi.resize(Dimension, nCorners);

    double csi = xi[0];

    phi[0] = (1. - csi)/2.;
    phi[1] = (1. + csi)/2.;

    dphi(0,0) = -1./2.;
    dphi(0,1) =  1./2.;

}

void Geom1d::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() < NodeCo.rows()) DebugStop(); // x.size() != NodeCo.rows()
    if(NodeCo.cols() !=nCorners) DebugStop();

    int nrow=NodeCo.rows();

    VecDouble phi(2);
    MatrixDouble dphi(Dimension,2);

    Shape(xi, phi, dphi);

    for (int i=0; i < nrow; i++) {
        x[i] = NodeCo(i,0)*phi[0] + NodeCo(i,1)*phi[1];
//        x[i] = NodeCo(i,0)*(1.-xi[0])*(1./2.) + NodeCo(i,1)*(1. + xi[0])*1./2.;
    }
 
}

void Geom1d::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() < NodeCo.rows()) DebugStop();  // x.size() != NodeCo.rows()
    if(NodeCo.cols() !=nCorners) DebugStop();

    int nrow=NodeCo.rows();
    int ncol=NodeCo.cols();

    gradx.resize(nrow, Dimension);
    gradx.setZero();
    x.resize(nrow);
    x.setZero();

    VecDouble phi(nCorners); // porque nCorners == 2
    MatrixDouble dphi(Dimension, nCorners);
    Shape(xi, phi, dphi);
//    X(xi, NodeCo, x);

    for(int i = 0; i < ncol; i++){
        for(int j = 0; j < nrow; j++){
            x[j] += NodeCo(j,i) * phi[i];
            gradx(j, 0) += NodeCo(j, i) * dphi(0, i);       
        } 
    }

}

void Geom1d::SetNodes(const VecInt &nodes) {
    if(nodes.rows() != 2) DebugStop();
    fNodeIndices = nodes;
}

void Geom1d::GetNodes(VecInt &nodes) const {
    nodes = fNodeIndices;
}


int Geom1d::NodeIndex(int node) const{
    // implementando a seguinte linha como extra
    if(node < 0 || node > 2) DebugStop();
    return fNodeIndices[node];
}

int Geom1d::NumNodes(){
    return nCorners;    
}


GeoElementSide Geom1d::Neighbour(int side) const {
    // implementando a seguinte linha como extra
    if(side < 0 || side > 2) DebugStop();
    return fNeighbours[side];
}

void Geom1d::SetNeighbour(int side, const GeoElementSide &neighbour) {
    // implementando a seguinte linha como extra
    if(side < 0 || side > 2) DebugStop();
    fNeighbours[side]=neighbour;
}
