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
#include "GeoMesh.h"
#include "ReadGmsh.h"
#include "CompMesh.h"
#include "VTKGeoMesh.h"
#include "Poisson.h"
#include "L2Projection.h"
#include "Analysis.h"
#include "PostProcessTemplate.h"

// verificar se os headers abaixo sao necessarios
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"

#include "CompElement.h"
#include "GeoElement.h"
#include "Assemble.h"
// ++++++++++++++++++++++++++++++++++++++++++++++

#ifdef WIN32
#define __PRETTY_FUNCTION__ __FUNCTION__
#endif

using std::cout;
using std::endl;
using std::cin;

int main ()
{
    GeoMesh gmesh;  // ler a malha que criamos
    ReadGmsh read;

    // std::string filename("quads1elements.msh");
    // std::string filename("quads4elements.msh");
    std::string filename("quads16elements.msh");
    // std::string filename("quads64elements.msh");    
    // std::string filename("quads16x16elements.msh");
    // std::string filename("quads32x32elements.msh");       

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
    cmesh.SetDefaultOrder(2); // ordem de aproximacao
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
    // postprocess.AppendVariable("DSol");
    postprocess.AppendVariable("SolExact");
    // postprocess.AppendVariable("Force");
    // postprocess.AppendVariable("DSolExact");    

    postprocess.SetExact(exact);
    mat1->SetExactSolution(exact);

    VTKGeoMesh plotmesh;
    // plotmesh.PrintGMeshVTK(&gmesh, "g_quads36elements.vtk"); 
    plotmesh.PrintGMeshVTK(&gmesh, "gquads4x4elementsorder2.vtk"); 
    // AnalysisLoc.PostProcessSolution("cquads36elements.vtk", postprocess);

    // Gerando a malha computacional para elementos quadrilateros de ordem 1
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "cquads1x1elements.vtk");    
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "cquads2x2elements.vtk");    
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "cquads4x4elements.vtk");                   
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "cquads8x8elements.vtk");
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "cquads16x16elements.vtk");
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "cquads32x32elements.vtk");
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "cquads64x64elements.vtk");


    // Gerando a malha computacional para elementos quadrilateros de ordem 2
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "cquads1x1elementsorder2.vtk");    
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "cquads2x2elementsorder2.vtk");      
    plotmesh.PrintCMeshVTK(&cmesh, 2, "cquads4x4elementsorder2.vtk");                 
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "cquads8x8elementsorder2.vtk");
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "cquads16x16elementsorder2.vtk");
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "cquads32x32elementsorder2.vtk");

    // Lembrete: neste codigo lemos o numero de elementos como:
    //           elements = side x side + 4 x side + 4 pts (nodes)
    // ----------------------------------------------------------

    VecDouble errvec;
    errvec = AnalysisLoc.PostProcessError(std::cout, postprocess); 
    
    return 0;
}

//     // std::vector<MathStatement *> mathvec = {0,mat1,bc_point,bc_linha};
//     // cmesh.SetMathVec(mathvec);
//     // cmesh.SetDefaultOrder(1);
//     // cmesh.AutoBuild();
//     // cmesh.Resequence();

//     //     Analysis locAnalysis(&cmesh);
//     // locAnalysis.RunSimulation();
//     // PostProcessTemplate<Poisson> postprocess;
//     // auto exact = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
//     // {
//     //     val[0] = (1.-x[0])*x[0]*(1-x[1])*x[1];
//     //     deriv(0,0) = (1.-2.*x[0])*(1-x[1])*x[1];
//     //     deriv(1,0) = (1-2.*x[1])*(1-x[0])*x[0];
//     // };

// //    if (!strcmp("Sol", name.c_str())) return ESol;
// //    if (!strcmp("DSol", name.c_str())) return EDSol;
// //    if (!strcmp("Flux", name.c_str())) return EFlux;
// //    if (!strcmp("Force", name.c_str())) return EForce;
// //    if (!strcmp("SolExact", name.c_str())) return ESolExact;
// //    if (!strcmp("DSolExact", name.c_str())) return EDSolExact;
//     postprocess.AppendVariable("Sol");
//     postprocess.AppendVariable("DSol");
//     postprocess.AppendVariable("Flux");
//     postprocess.AppendVariable("Force");
//     postprocess.AppendVariable("SolExact");
//     postprocess.AppendVariable("DSolExact");
//     postprocess.SetExact(exact);
//     mat1->SetExactSolution(exact);
//     locAnalysis.PostProcessSolution("quads.vtk", postprocess);

//     VecDouble errvec;
//     errvec = locAnalysis.PostProcessError(std::cout, postprocess);
    
//     return 0;
// }