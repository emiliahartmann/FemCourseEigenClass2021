

#include <iostream>
#include <math.h>
#include "IntRule.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTetrahedron.h"
#include "IntRuleTriangle.h"
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "DataTypes.h"
#include "Analysis.h"
#include "VTKGeoMesh.h"
#include "ReadGmsh.h"
#include "CompMesh.h"
#include "Poisson.h"
#include "L2Projection.h"   // vai ser usado para função de contorno
//#include "CompElement.h"
//#include "MathStatement.h"

using std::cout;
using std::endl;
using std::cin;

int main (){
    VecDouble phir(2),phitheta(2);
    MatrixDouble dphir(1,2),dphitheta(1,2);
    VecDouble xp;
    int order = 1;
    MatrixDouble jac(2,2);

    double Integral = 0;
    IntRuleQuad rule(1);
    int np = rule.NPoints();
    for(int ip=0; ip<np; ip++)
    {
        VecDouble xip(2);
        double wp;
        rule.Point(ip, xip, wp);
        
        double r = (xip[0]-1.)/2.+5.*(xip[0]+1.)/2.;
        double drdxi = 1./2.+5./2.;
        double theta = M_PI/2.*(xip[1]+1.)/2.;
        double dthetadeta = M_PI/4.;
        // x = r cos(theta)
        // y = r sin(theta)
        // dxdr = r' cos(theta)
        jac(0,0) = drdxi*cos(theta);
        // dxdtheta = -r sin(theta) theta'
        jac(0,1) = -r*sin(theta)*dthetadeta;
        // dydr = r' sin(theta)
        jac(1,0) = drdxi*sin(theta);
        // dydtheta = r cos(theta) theta'
        jac(1,1) = r*cos(theta)*dthetadeta;
        double detjac = std::abs(jac.determinant());
        Integral += detjac * wp;
    }
    
    std::cout << "order = " << order << " integral aproximada " << Integral <<
    " erro " << 6.*M_PI-Integral << std::endl;


    GeoMesh gmesh;          // ler a malha que criamos
        ReadGmsh read;
        read.Read(gmesh,"/home/emilia/Repositórios/FemCourseEigenClass2021/mainprograms/malha.msh");
        VTKGeoMesh plotmesh;
        plotmesh.PrintGMeshVTK(&gmesh, "malha.vtk");

    CompMesh cmesh(&gmesh);  // criar uma malha computacional
        MatrixDouble perm(2,2);
        perm.setZero();
        perm(0,0) = 1.;
        perm(1,1) = 1.;
        Poisson *mat1 = new Poisson(1, perm);

        MatrixDouble proj(1,1), val1(1,1), val2(1,1);
        proj.setZero();
        val1.setZero();
        // val2.setZero();
        val2.setOnes();
        L2Projection *bc_linha = new L2Projection(0, 2, proj, val1, val2); // 0 -> 1 para condicao de Neumann
        L2Projection *bc_point = new L2Projection(0, 3, proj, val1, val2);
        std::vector<MathStatement*>mathvec (4);
        mathvec[0] = 0;
        mathvec[1] = mat1;
        mathvec[2] = bc_linha;
        mathvec[3] = bc_point;
        //= {0, mat1, bc_linha, bc_point}; // o meu material tem id igual a 1
        cmesh.SetMathVec(mathvec); 
        cmesh.AutoBuild();
        cmesh.Resequence();
        cmesh.Solution() (100,0) = 1.;    // apenas para teste, excluir esta linha depois -- para visualizar as funcoes de forma
        plotmesh.PrintCMeshVTK(&cmesh, 2, "cmesh_malha.vtk"); 
    return 0;

}
