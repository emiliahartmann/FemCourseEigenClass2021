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
#include "CompMesh.h"
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

void CreateTestMesh(CompMesh &mesh, int order, double h);

void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv);

int main ()
{
    
    
    CompMesh mesh;
    int order = 3;
    double h=1.0/8.0;
    CreateTestMesh(mesh,order,h);
    
    Analysis Analysis(&mesh);
    Analysis.RunSimulation();
    
    //exact(x, val, deriv);
    
    double energy=0.0,l2=0.0;
    
    PostProcessTemplate<Poisson> postprocess;
    postprocess.SetExact(exact);
    
    VecDouble errvec;
    errvec = Analysis.PostProcessError(std::cout, postprocess);
    
    cout<<endl;
    cout<<energy<<endl;
    cout<<l2<<endl;
    

    //TestOneDProblem(&mesh);
    
    return 0;
}
void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv){


    double E=exp(1.0);
    VecDouble x(1);
    x[0]=point[0];
    
    val[0]=(30. + 100.*pow(E,100.) - 130.*pow(E,10.*x[0]) - 3*x[0] + 3*pow(E,100.)*x[0])/(10.*(-1. + pow(E,100.)));
    deriv(0,0)=(-3. + 3*pow(E,100) - 1300*pow(E,10*x[0]))/(10.*(-1 + pow(E,100)));
}


void CreateTestMesh(CompMesh &mesh, int order, double h)
{
    DebugStop();
}