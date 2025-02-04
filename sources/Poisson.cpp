/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Poisson.h"
///\cond
#include <functional>
#include <string.h>
///\endcond

Poisson::Poisson() {
}

Poisson::Poisson(int materialid, MatrixDouble &perm) {
    permeability = perm;
    this->SetMatID(materialid);
}

Poisson::Poisson(const Poisson &copy) {
    permeability = copy.permeability;
    forceFunction = copy.forceFunction;
}

Poisson &Poisson::operator=(const Poisson &copy) {
    permeability = copy.permeability;
    forceFunction = copy.forceFunction;
    return *this;
}

Poisson *Poisson::Clone() const {
    return new Poisson(*this);
}

Poisson::~Poisson() {
}

MatrixDouble Poisson::GetPermeability() const {
    return permeability;
}

void Poisson::SetPermeability(const MatrixDouble &perm) {
    permeability = perm;
}

int Poisson::NEvalErrors() const {
    return 3;
}

int Poisson::VariableIndex(const PostProcVar var) const {
    return int(var);
}

Poisson::PostProcVar Poisson::VariableIndex(const std::string &name) {
    if (!strcmp("Sol", name.c_str())) return ESol;
    if (!strcmp("DSol", name.c_str())) return EDSol;
    if (!strcmp("Flux", name.c_str())) return EFlux;
    if (!strcmp("Force", name.c_str())) return EForce;
    if (!strcmp("SolExact", name.c_str())) return ESolExact;
    if (!strcmp("DSolExact", name.c_str())) return EDSolExact;
    else {
        std::cout << "variable "<< name << " not implemented" << std::endl;
    }
    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return ENone;
}

int Poisson::NSolutionVariables(const PostProcVar var) {
    if (var == ESol) return this->NState();
    if (var == EDSol) return 3;
    if (var == EFlux) return 3;
    if (var == EForce) return this->NState();
    if (var == ESolExact) return this->NState();
    if (var == EDSolExact) return 3;
    else {
        std::cout << "variable not implemented" << std::endl;
    }
    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return -1;
}

void Poisson::ContributeError(IntPointData &data, VecDouble &u_exact, MatrixDouble &du_exact, VecDouble &errors) const {
    errors.resize(NEvalErrors());
    errors.setZero();
    MatrixDouble gradu;
    MatrixDouble axes = data.axes;

    VecDouble u = data.solution;
    MatrixDouble dudx = data.dsoldx; // futuramente, mudar o nome dsoldx para dsoldaxes

    this->Axes2XYZ(dudx, gradu, axes);

    double diff = 0.0;
    for (int i = 0; i < this->NState(); i++) {
        diff = (u[i] - u_exact[i]);
        errors[0] += diff*diff;
    } // Erro L2

    errors[1] = 0.;
    int dim = this->Dimension();
    int nstate = this->NState();
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < nstate; j++) {
            diff = (gradu(i, j) - du_exact(i, j));
            errors[1] += diff*diff;
        } // aqui o erro deveria ser diff*perm*diff..Num codigo com diferentes permeabilidades, teriamos um problema

    } // Erro da derivada: h1 c menor
    errors[2] = errors[0] + errors[1]; // norma do erro: produto interno da funcao e de sua derivada
}

void Poisson::Contribute(IntPointData &data, double weight, MatrixDouble &EK, MatrixDouble &EF) const {

    VecDouble phi = data.phi;
    MatrixDouble dphi = data.dphidx; // dphidx = dphidaxis na vdd
    MatrixDouble axes = data.axes;   // por isso que precisamos reescrever os dphidx por dphi2 e dphi3
    MatrixDouble dphi2, dphi3;

    // versao do professor
    int nstate = this->NState();
    dphi2 = data.axes.transpose()*data.dphidx;
    dphi3 = dphi2.transpose();

    MatrixDouble perm(3, 3);
    perm = this->GetPermeability();
    double res = 0.; // vamos fazer uma aproximacao aqui com menos laplaciano de u = fx (eq de Poisson)
                     // isso pode ser feito colocando um ponteiro que calcula fx 

    auto force = this->GetForceFunction();
    if(force){
        VecDouble resloc(nstate);
        force(data.x, resloc);
        res = resloc[0]; 
    }

    // res = 1.; // Temporario, apenas para visualizar o vetor de carga diferente de zero

    // std::cout << "perm/n" << perm << std::endl;
    // std::cout << "dphi2/n" << dphi2 << std::endl;    
    // std::cout << "dphi3/n" << dphi3 << std::endl;
    // std::cout << "res/n" << res << std::endl;    

    // eq. diferencial = -u'' = fx
    // EF += phi*(res*weight);
    // EK += dphi3*perm*dphi2*weight;

    // eq. diferencial = -u''+ u = fx
    EF += phi*(res*weight);
    EK += (dphi3*perm*dphi2 + phi*phi.transpose())*weight;

    // versao a Emilia
    // this->Axes2XYZ(dphi, dphi2, axes);
    // dphi3 = dphi2.transpose();

    // int nshape = phi.size();
    // int nstate = this->NState();
    // int dim = dphi.rows();

    // MatrixDouble perm(3, 3);
    // std::function<void(const VecDouble &co, VecDouble & result) > force;

    // perm = this->GetPermeability();
    // double resloc = 0.;
    // force = this->GetForceFunction();

    
    // if(force)
    // {
    //     VecDouble res(nstate);
    //     force(data.x, res);
    //     resloc = res[0];
    // }

    // //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // // Fazer o calculo da matriz de rigidez e da carga de forca aqui 
    // for (int i = 0; i < nshape; i++) {
    //     for (int j = 0; j < nshape; j++) {
    //         for (int d=0; d < dim; d++) {
    //             EK(i,j) += dphi2(d,i)*dphi3(d,j)*weight*perm(i,j);
    //         }
    //         EF(i,0) += phi(i,0)*weight*resloc;
    //     // std::cout << "EK( " << i << "," << j << " ) " << EK(i,j) << " " << std::endl;  
    //     }
        
    // } 

    //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}

void Poisson::PostProcessSolution(const IntPointData &data, const int var, VecDouble &Solout) const {
    // VecDouble sol = data.solution;
    // int solsize = sol.size();
    // int rows = data.dsoldx.rows();  // nao precisa mais
    // int cols = data.dsoldx.cols();  // nao precisa mais
    // MatrixDouble gradu(rows, cols); // nao precisa mais
    // gradu = data.dsoldx;            // nao precisa mais  

    MatrixDouble gradudx, flux;
    gradudx = data.axes.transpose()*data.dsoldx;
    flux = -permeability*gradudx;

    int nstate = this->NState();
    if(nstate != 1) DebugStop();

    switch (var) {
        case 0: //None
        {
            std::cout << " Var index not implemented " << std::endl;
            DebugStop();

        }

        case 1: //ESol
        {
            Solout = data.solution;
        }
            break;

        case 2: //EDSol
        {
            //+++++++++++++++++
            // versao local
            // Solout.resize(rows, cols);
            // Solout = gradu;

            // versao do professor
            Solout.resize(3);
            for (int i = 0; i < 3; i++){
                Solout[i] = gradudx(i, 0);
            }
            //+++++++++++++++++           
        }
            break;
        case 3: //EFlux
        {
            //+++++++++++++++++
            Solout.resize(3);
            for (int i = 0; i < 3; i++){
                Solout[i] = flux(i, 0);
            }
            //+++++++++++++++++
        }
            break;

        case 4: //EForce
        {   // versao local
            // Solout.resize(nstate);
            // VecDouble result(nstate);
            // this->forceFunction(data.x, result);
            // for (int i = 0; i < nstate; i++) {
            //     Solout[i] = result[i];
            // }

            // versao do professor
            Solout.resize(nstate);
            if(forceFunction) this->forceFunction(data.x, Solout);
            else Solout.setZero();
        }
            break;

        case 5: //ESolExact
        {
            //+++++++++++++++++
            Solout.resize(nstate);
            VecDouble sol(nstate);
            MatrixDouble dsol(3, nstate);
            if(SolutionExact) this->SolutionExact(data.x, Solout, dsol);
            else Solout.setZero();
            //+++++++++++++++++
        }
            break;
        case 6: //EDSolExact
        {
            //+++++++++++++++++
            Solout.resize(3);
            VecDouble sol(nstate);
            MatrixDouble dsol(3, nstate);
            if(SolutionExact) this->SolutionExact(data.x, sol, dsol);
            else dsol.setZero();

            for (int i = 0; i < 3; i++) {
                Solout[i] = dsol(i, 0);
            }
            //+++++++++++++++++
        }
            break;
        default:
        {
            std::cout << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }   
}

double Poisson::Inner(MatrixDouble &S, MatrixDouble & T) const {
    double inner = 0;
    for (int i = 0; i < S.rows(); i++) {
        for (int j = 0; j < S.cols(); j++) {
            inner += S(i, j) * T(i, j);
        }
    }
    return inner;
}
