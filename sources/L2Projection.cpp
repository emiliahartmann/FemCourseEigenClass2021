/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "L2Projection.h"
#include "PostProcess.h"
#include "tpanic.h"
#include <string.h>

L2Projection::L2Projection() {
}

L2Projection::L2Projection(int bctype, int materialid, MatrixDouble &proj, MatrixDouble Val1, MatrixDouble Val2) {
    projection = proj;
    BCType = bctype;
    BCVal1 = Val1;
    BCVal2 = Val2;
    this->SetMatID(materialid);
}

L2Projection::L2Projection(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
    BCType = copy.BCType;
    BCVal1 = copy.BCVal1;
    BCVal2 = copy.BCVal2;

}

L2Projection &L2Projection::operator=(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
    BCType = copy.BCType;
    BCVal1 = copy.BCVal1;
    BCVal2 = copy.BCVal2;
    return *this;
}

L2Projection *L2Projection::Clone() const {
    return new L2Projection(*this);
}

L2Projection::~L2Projection() {
}

MatrixDouble L2Projection::GetProjectionMatrix() const {
    return projection;
}

void L2Projection::SetProjectionMatrix(const MatrixDouble &proj) {
    projection = proj;
}

void L2Projection::Contribute(IntPointData &data, double weight, MatrixDouble &EK, MatrixDouble &EF) const {
    int nstate = this->NState();
    int nshape = data.phi.size();

    VecDouble result(data.x.size());
    MatrixDouble deriv(data.x.size(), data.x.size());

    SolutionExact(data.x, result, deriv);

    switch (this->GetBCType()) {

        case 0:
        {
            for (int iv = 0; iv < nstate; iv++) {
                for (int in = 0; in < nshape; in++) {
                    EF(nstate * in + iv, 0) += MathStatement::gBigNumber * result[iv] * data.phi[in] * weight;
                    for (int jn = 0; jn < nshape; jn++) {
                        EK(nstate * in + iv, nstate * jn + iv) += MathStatement::gBigNumber * data.phi[in] * data.phi[jn] * weight;
                    }
                }
            }
            break;
        }

        case 1:
        {
            for (int iv = 0; iv < nstate; iv++) {
                for (int in = 0; in < nshape; in++) {
                    EF(nstate * in + iv, 0) += Val2()(iv, 0) * data.phi[in] * weight;
                }
            }
            break;
        }

        default:
        {
            std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
        }
    }
}

int L2Projection::NEvalErrors() const {
    return 3;
}

void L2Projection::ContributeError(IntPointData &data, VecDouble &errors) const {
    return;
}

int L2Projection::VariableIndex(const PostProcVar var) const {
    if (var == ESol) return ESol;

    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return 0;
}

L2Projection::PostProcVar L2Projection::VariableIndex(const std::string & name) {
    if (!strcmp("Solution", name.c_str())) return ESol;

    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return ENone;
}

int L2Projection::NSolutionVariables(const PostProcVar var) {
    if (var == ESol) return this->NState();

    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return 0;
}

void L2Projection::PostProcessSolution(const IntPointData &integrationpointdata, const int var, VecDouble &sol) const {
  sol.resize(this->NState(), 0);
}
