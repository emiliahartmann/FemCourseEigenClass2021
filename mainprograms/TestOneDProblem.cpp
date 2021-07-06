//
//  TestOneDProblem.cpp
//  FemSC
//
//  Created by Eduardo Ferri on 08/17/15.
//
//
//TestOneDProblem cpp
/*
        Os testes foram preparados com um proposito educacional,
        recomenda-se que o aluno entenda a funcionalidade de cada
        teste e posteriormente use com seu código caso a caso
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
#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "MathStatement.h"
#include "L2Projection.h"
#include "Analysis.h"
#include "IntRule.h"
#include "PostProcessTemplate.h"
#include "Poisson.h"

using std::cout;
using std::endl;
using std::cin;

void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv); // declaracao de uma funcao da antiga linguagem C

int main ()
{
    GeoMesh gmesh;
    ReadGmsh read;
    // std::string filename("oned4elements.msh"); // 4 elements    
    // std::string filename("oned8elements.msh"); // 8 elements
    // std::string filename("oned16elements.msh"); // 16 elements
    // std::string filename("oned32elements.msh"); // 32 elements
    // std::string filename("oned80elements.msh"); // 80 elements   
    std::string filename("oned160elements.msh"); // 160 elements  
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
    mat1->SetDimension(1);

    auto force = [](const VecDouble &loc, VecDouble &f){
        const auto &x = loc[0];
        const auto &y = loc[1];

        f[0] = x;
        // f[0] = 1.;
    };   
            

    // auto force = [](const VecDouble &x, VecDouble &res)
    // {
    //     // res[0] = 1.; //x[0]
    //     res[0] = x[0]; 
    // }; // declaracao de funcao tbm; dessa forma temos mais autonomia na funcao que escolhemos.
        // os colchetes servem para setar quais variaveis do codigo vc vai usar para para construir a funcao force
    mat1->SetForceFunction(force);
    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    L2Projection *bc_linha = new L2Projection(0,2,proj,val1,val2);
    L2Projection *bc_point = new L2Projection(0,3,proj,val1,val2);
    std::vector<MathStatement *> mathvec = {0,mat1,bc_point,bc_linha};
    cmesh.SetMathVec(mathvec);
    cmesh.SetDefaultOrder(1); // ordem de aproximacao
    cmesh.AutoBuild();
    cmesh.Resequence();

    Analysis AnalysisLoc(&cmesh);
    AnalysisLoc.RunSimulation();
    
    PostProcessTemplate<Poisson> postprocess;
    postprocess.SetExact(exact);
    
    VecDouble errvec;
    errvec = AnalysisLoc.PostProcessError(std::cout, postprocess); 
    // Tomar cuidado pq no PostProcessaento, ele calcula o erro para o mathId igual a 1 apenas....tentar modificar isso futuramente
    
    // Imprimir a malha VTK
    VTKGeoMesh plotmesh;
    // Para ordem de aproximacao igual a 1
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "oned4elementsorder1.vtk");  
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "oned8elementsorder1.vtk");   
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "oned16elementsorder1.vtk");
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "oned32elementsorder1.vtk");   
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "oned80elementsorder1.vtk");
    plotmesh.PrintCMeshVTK(&cmesh, 2, "oned160elementsorder1.vtk");
    
    // Para ordem de aproximacao igual a 2 
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "oned4elements.vtk");   
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "oned8elements.vtk");     
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "oned16elements.vtk");    
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "oned32elements.vtk");   
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "oned80elements.vtk");
    // plotmesh.PrintCMeshVTK(&cmesh, 2, "oned160elements.vtk");
    return 0;
}
void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv){

    //// Para Laplaciano(u) = 1.;
    // deriv(0,0) = 4-point[0];
    // val[0]=point[0]*(8.-point[0])/2.;

    // // Para Laplaciano(u) = x;
    // deriv(0,0) = 32./3.-(point[0]*point[0]/2.);
    // val[0]=point[0]*(32.-point[0]*point[0]/2.)/3.;

    // Para u" + u = x;
    deriv(0,0) = 1. - (8./((exp(8.)-exp(-8.))))*(exp(point[0])+exp(-point[0]));
    val[0]= -(8./((exp(8.)-exp(-8.))))*(exp(point[0])-exp(-point[0]))+point[0];
    // double E=exp(1.0);
    // VecDouble x(1);
    // x[0]=point[0];
    
    // val[0]=(30. + 100.*pow(E,100.) - 130.*pow(E,10.*x[0]) - 3*x[0] + 3*pow(E,100.)*x[0])/(10.*(-1. + pow(E,100.)));
    // deriv(0,0)=(-3. + 3*pow(E,100) - 1300*pow(E,10*x[0]))/(10.*(-1 + pow(E,100)));



    return;
}


