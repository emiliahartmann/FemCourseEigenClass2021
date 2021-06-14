

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
        val2.setZero();
        L2Projection *bc_linha = new L2Projection(0, 2, proj, val1, val2);
        L2Projection *bc_point = new L2Projection(0, 11, proj, val1, val2);
        std::vector<MathStatement*>mathvec (12);
        mathvec[1] = mat1;
        mathvec[2] = bc_linha;
        mathvec[11] = bc_point;
        //= {0, mat1, bc_linha, bc_point}; // o meu material tem id igual a 1
        cmesh.SetMathVec(mathvec); 
        cmesh.AutoBuild();
        cmesh.Resequence();
        cmesh.Solution() (0,0) = 1.;
        plotmesh.PrintCMeshVTK(&cmesh, 2, "cmesh_malha.vtk");
    return 0;

}
