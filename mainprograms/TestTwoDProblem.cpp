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
//#include "GeoNode.h"
//#include "GeoElement.h"
//#include "IntPointData.h"
#include "CompElementTemplate.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "CompMesh.h"
#include "GeoMesh.h"
//#include "GeoElement.h"
//#include "GeoElementTemplate.h"
//#include "MathStatement.h"
#include "Poisson.h"
//#include "L2Projection.h"
//#include "Analysis.h"
//#include "IntRule.h"
//#include "PostProcess.h"
#include "PostProcessTemplate.h"

#ifdef WIN32
#define __PRETTY_FUNCTION__ __FUNCTION__
#endif

class CompMesh;

using std::cout;
using std::endl;
using std::cin;

void CreateTestMesh(CompMesh &mesh, int order);

void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv);

int main ()
{
    
    
    CompMesh mesh;
    int order = 2;
    //double h=1.; //para dar el tamanho del elemento
    CreateTestMesh(mesh,order);

    mesh.Print();

	Analysis Analisis(&mesh);
	Analisis.RunSimulation();

    VecDouble vecerr;
    PostProcessTemplate<Poisson> postprocess;
    postprocess.SetExact(exact);
    vecerr = Analisis.PostProcessError(std::cout, postprocess);
    
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


void CreateTestMesh(CompMesh &mesh, int order)
{
    DebugStop();
}