//
//  TestOneDProblem.cpp MODIFICADO DO ORIGINAL
//  FemSC
//
//  Created by Eduardo Ferri on 08/17/15.
//
//
//TestOneDProblem cpp
/*
 Os testes foram preparados com um proposito educacional,
 recomenda-se que o aluno entenda a funcionalidade de cada
 teste e posteriormente use com seu cÛdigo caso a caso
 */
//      Obs: O xmax e xmin estao tomados como 4 e 0, respectivamente,
//      caso estes valores sejam alterados, editar o teste TestNodes.
//
//
#include <iostream>
#include <math.h>

#include "CompElementTemplate.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "CompMesh.h"
#include "GeoMesh.h"
#include "ReadGmsh.h"
#include "VTKGeoMesh.h"
#include "Poisson.h"
#include "L2Projection.h"
#include "Analysis.h"
#include "PostProcessTemplate.h"

#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"

#include "CompElement.h"
#include "GeoElement.h"
#include "Assemble.h"

#ifdef WIN32
#define __PRETTY_FUNCTION__ __FUNCTION__
#endif

// class CompMesh;

using std::cout;
using std::endl;
using std::cin;

// void CreateTestMesh(CompMesh &mesh, int order);

// void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv);



int main ()
{
    GeoMesh gmesh;  // ler a malha que criamos
    ReadGmsh read;
    std::string filename("quads36elements.msh");
#ifdef MACOSX
    filename = "../"+filename;
#endif
    read.Read(gmesh,filename);
    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) = 1.;
    perm(1,1) = 1.;
    perm(2,2) = 1.;
    Poisson *mat1 = new Poisson(1,perm);
    mat1->SetDimension(2); //2
    
    auto force = [](const VecDouble &loc, VecDouble &f)
{
    const auto &x = loc[0];
    const auto &y = loc[1];

    // f[0] = 2.*(1.-x)*x + 2.*(1.-y)*y; // para o problema -lap(u) = fx
    f[0] = +2.*(1.-x)*x + 2.*(1.-y)*y + (1.-x)*x*(1.-y)*y; // para o problema -lap(u) +u = fx
    
    // res[0] = 2.*(1-x[0])*x[0] + 2.*(1-x[1])*x[1]; 
};

    mat1->SetForceFunction(force);
    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    L2Projection *bc_linha = new L2Projection(0,2,proj,val1,val2);
    L2Projection *bc_point = new L2Projection(0,3,proj,val1,val2);
    // bc_linha->SetExactSolution(exact);
    // bc_point->SetExactSolution(exact);
    std::vector<MathStatement *> mathvec = {0,mat1,bc_linha,bc_point};
    cmesh.SetMathVec(mathvec);
    cmesh.SetDefaultOrder(1); // ordem de aproximacao
    cmesh.AutoBuild();
    cmesh.Resequence();

    Analysis AnalysisLoc(&cmesh);
    AnalysisLoc.RunSimulation();
    
    PostProcessTemplate<Poisson> postprocess;
    
    auto exact = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        // val[0] = (1-x[0])*x[0]+(1-x[1])*x[1]; 
        // deriv(0,0) = (1-2.*x[0]); 
        // deriv(1,0) = (1-2.*x[1]);         
        val[0] = (1-x[0])*x[0]*(1-x[1])*x[1]; 
        deriv(0,0) = (1-2.*x[0])*(1-x[1])*x[1]; 
        deriv(1,0) = (1-x[0])*x[0]*(1-2.*x[1]); 
    };

    // postprocess.AppendVariable("Flux");
    // postprocess.AppendVariable("Sol");
    postprocess.AppendVariable("DSol");
    // postprocess.AppendVariable("SolExact");
    // postprocess.AppendVariable("Force");
    // postprocess.AppendVariable("DSolExact");    

    postprocess.SetExact(exact);
    mat1->SetExactSolution(exact);

    VTKGeoMesh plotmesh;
    plotmesh.PrintGMeshVTK(&gmesh, "g_quads36elements.vtk"); 
    plotmesh.PrintCMeshVTK(&cmesh, 2, "c_quads36elements.vtk"); 
    // AnalysisLoc.PostProcessSolution("quads36elements.vtk", postprocess);

    VecDouble errvec;
    errvec = AnalysisLoc.PostProcessError(std::cout, postprocess); 
    
    return 0;
}


void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv){
    
    
    long double Pi=4*atan(1);
    int mmax=200;
    int nmax=200;
    
//    val=0.;
//    deriv[0]=0.;
//    deriv[1]=0.;
//    
//    val=-point[0]*point[0]+point[0];
//    deriv[0]=-2.*point[0]+1.;
//    deriv[1]=0.;
    
    for(int m=1;m<=mmax;m++){
        for(int n=1;n<=nmax;n++){
            val[0]+=(4*(-1 + pow(-1,m))*(-1 + pow(-1,n))*sin(m*Pi*point[0])*sin(n*Pi*point[1]))/(m*n*(pow(m,2) + pow(n,2))*pow(Pi,4));
            
            deriv(0,0)+=(4*(-1 + pow(-1,m))*(-1 + pow(-1,n))*cos(m*Pi*point[0])*sin(n*Pi*point[1]))/(n*(pow(m,2) + pow(n,2))*pow(Pi,3));
            
            deriv(1,0)+=(4*(-1 + pow(-1,m))*(-1 + pow(-1,n))*cos(n*Pi*point[1])*sin(m*Pi*point[0]))/(m*(pow(m,2) + pow(n,2))*pow(Pi,3));
            
        }
    }
    
    
    
    
}


// void CreateTestMesh(CompMesh &mesh, int order)
// {
//     DebugStop();
// }