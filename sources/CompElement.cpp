/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "CompElement.h"
#include "GeoElement.h"
#include "MathStatement.h"
#include "CompMesh.h"
#include "tpanic.h"
#include <math.h> 
using namespace std;

CompElement::CompElement() {
}

CompElement::CompElement(int64_t ind, CompMesh *cmesh, GeoElement *geo) {
    compmesh = cmesh;
    index = ind;
    geoel = geo;
}

CompElement::CompElement(const CompElement &copy) {
    compmesh = copy.compmesh;
    index = copy.index;
    geoel = copy.geoel;
    intrule = copy.intrule;
    mat = copy.mat;
}

CompElement &CompElement::operator=(const CompElement &copy) {
    compmesh = copy.compmesh;
    index = copy.index;
    geoel = copy.geoel;
    intrule = copy.intrule;
    mat = copy.mat;
    return *this;
}

CompElement::~CompElement() {
}

// CompElement *CompElement::Clone() const {
//     return new CompElement(*this);
// }

DOF &CompElement::GetDOF(int i) const
{
    return compmesh->GetDOF(GetDOFIndex(i));
}

MathStatement *CompElement::GetStatement() const {
    return mat;
}

void CompElement::SetStatement(MathStatement *statement) {
    mat = statement;
}

IntRule *CompElement::GetIntRule() const {
    return intrule;
}

void CompElement::SetIntRule(IntRule *irule) {
    intrule = irule;
}

void CompElement::SetIndex(int64_t ind) {
    index = ind;

}

GeoElement *CompElement::GetGeoElement() const {
    return geoel;
}

void CompElement::SetGeoElement(GeoElement *element) {
    geoel = element;
}

CompMesh *CompElement::GetCompMesh() const {
    return compmesh;
}

void CompElement::SetCompMesh(CompMesh *mesh) {
    compmesh = mesh;
}

void CompElement::InitializeIntPointData(IntPointData &data) const {

    const int dim = this->Dimension();
    const int nshape = this->NShapeFunctions();
    const int nstate = this->GetStatement()->NState();

    data.weight = 0;
    data.detjac = 0;
    data.phi.resize(nshape, 0);
    data.dphidx.resize(dim, nshape);
    data.dphidksi.resize(dim, nshape);
    data.x.resize(3);
    data.ksi.resize(dim, 0);
    data.gradx.resize(dim, nshape);
    data.axes.resize(dim, 3);
    data.solution.resize(1, 0);
    data.dsoldksi.resize(1, 1);
    data.dsoldx.resize(1, 1);

}

void CompElement::ComputeRequiredData(IntPointData &data, VecDouble &intpoint) const {
    data.ksi = intpoint;

    int dim = this->Dimension();

    geoel->X(intpoint, data.x);
    geoel->GradX(intpoint, data.x, data.gradx);

    data.detjac = 0.0;
    MatrixDouble jac(dim, dim);
    jac.setZero();
    MatrixDouble jacinv(dim, dim);
    jacinv.setZero();

    geoel->Jacobian(data.gradx, jac, data.axes, data.detjac, jacinv);

    this->ShapeFunctions(intpoint, data.phi, data.dphidksi);
    this->Convert2Axes(data.dphidksi, jacinv, data.dphidx);
}

void CompElement::Convert2Axes(const MatrixDouble &dphi, const MatrixDouble &jacinv, MatrixDouble &dphidx) const {
    int nshape = this->NShapeFunctions();
    int dim = this->Dimension();

    int ieq;
    switch (dim) {
        case 0:
        {

        }
            break;
        case 1:
        {
            for (ieq = 0; ieq < nshape; ieq++) {
                dphidx(0, ieq) *= jacinv(0, 0);
            }
        }
            break;
        case 2:
        {
            for (ieq = 0; ieq < nshape; ieq++) {
                dphidx(0, ieq) = jacinv(0, 0) * dphi(0, ieq) + jacinv(1, 0) * dphi(1, ieq);
                dphidx(1, ieq) = jacinv(0, 1) * dphi(0, ieq) + jacinv(1, 1) * dphi(1, ieq);
            }
        }
            break;
        case 3:
        {
            for (ieq = 0; ieq < nshape; ieq++) {
                dphidx(0, ieq) = jacinv(0, 0) * dphi(0, ieq) + jacinv(1, 0) * dphi(1, ieq) + jacinv(2, 0) * dphi(2, ieq);
                dphidx(1, ieq) = jacinv(0, 1) * dphi(0, ieq) + jacinv(1, 1) * dphi(1, ieq) + jacinv(2, 1) * dphi(2, ieq);
                dphidx(2, ieq) = jacinv(0, 2) * dphi(0, ieq) + jacinv(1, 2) * dphi(1, ieq) + jacinv(2, 2) * dphi(2, ieq);
            }
        }
            break;
    }
}

void CompElement::CalcStiff(MatrixDouble &ek, MatrixDouble &ef) const {
    MathStatement *material = this->GetStatement();
    if (!material) {
        std::cout << "Error at CompElement::CalcStiff" << std::endl;
        return;
    }
    
    ek.setZero();
    ef.setZero();

    IntPointData data;
    this->InitializeIntPointData(data);
    double weight = 0.;

    IntRule *intrule = this->GetIntRule();
    int intrulepoints = intrule->NPoints();    

    for (int int_ind = 0; int_ind < intrulepoints; ++int_ind) {
        std::cout << "Please insert the correct code\n";
        DebugStop();
    }
}

void CompElement::EvaluateError(VecDouble &errors) const {
    MathStatement * material = this->GetStatement();

    if (!material) {
        std::cout << "No material for this element\n";
        return;
    }

    int NErrors = material->NEvalErrors();
    errors.resize(NErrors, 0);


    IntRule *intrule = this->GetIntRule();
    int maxIntOrder = 15;
    intrule->SetOrder(maxIntOrder);

    int dim = Dimension();
    int nstate = material->NState();

    VecDouble values(NErrors, 1);
    double weight = 0.;

    IntPointData data;
    this->InitializeIntPointData(data);
    int nintpoints = intrule->NPoints();

    VecDouble u_exact(nstate);
    MatrixDouble du_exact(dim, nstate);

    for (int nint = 0; nint < nintpoints; nint++) {
        values.setZero();
        intrule->Point(nint, data.ksi, weight);
        this->ComputeRequiredData(data, data.ksi);
        weight *= fabs(data.detjac);

        this->GetMultiplyingCoeficients(data.coefs);
        data.ComputeSolution();
        material->ContributeError(data, values);

        for (int ier = 0; ier < NErrors; ier++)
            errors[ier] += weight * values[ier];
    }

    for (int ier = 0; ier < NErrors; ier++) {
        errors[ier] = sqrt(errors[ier]);
    }
}

void CompElement::Solution(VecDouble &intpoint, int var, VecDouble &sol) const {
    MathStatement * material = this->GetStatement();
    if (!material) {
        std::cout << "No material for this element\n";
        return;
    }

    IntPointData data;
    this->InitializeIntPointData(data);
    this->ComputeRequiredData(data, intpoint);
    this->GetMultiplyingCoeficients(data.coefs);
    data.ComputeSolution();

    material->PostProcessSolution(data, var, sol);
}
