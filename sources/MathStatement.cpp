/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "MathStatement.h"

double MathStatement::gBigNumber = 1.e12;

MathStatement::MathStatement() {
}

MathStatement::MathStatement(const MathStatement &copy) {
    matid = copy.matid;
}

MathStatement &MathStatement::operator=(const MathStatement &copy) {
    matid = copy.matid;
    return *this;
}

MathStatement::~MathStatement() {
}

void MathStatement::Axes2XYZ(const MatrixDouble &dudaxes, MatrixDouble &dudx, const MatrixDouble &axesv, bool colMajor) const {
    MatrixDouble axes(axesv.rows(), axesv.cols());
    for (int r = 0; r < axes.rows(); r++) {
        for (int c = 0; c < axes.cols(); c++) {
            axes(r, c) = axesv(r, c);
        }
    }

    if (colMajor) {
        MatrixDouble axesT;
        axesT = axes.transpose();

        if (dudx.rows() != axesT.rows() || dudx.cols() != dudaxes.cols()) {
            dudx.resize(axesT.rows(), dudaxes.cols());
        }
        dudx.setZero();
        dudx = axesT*dudaxes;
//        axesT.Multiply(dudaxes, dudx, 0);
    } else {
        dudx.resize(dudaxes.rows(), axes.cols());
        dudx.setZero();
        dudx = dudaxes*axes;
//        dudaxes.Multiply(axes, dudx, 0);
    }
}

void MathStatement::Print(std::ostream &out) {
    out << "Material: " << this->GetMatID() << std::endl;
    out << "Big number: " << gBigNumber << std::endl << std::endl;

}
