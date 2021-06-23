
///\cond
#include <iostream>
#include <math.h>
///\endcond
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "Analysis.h"
#include "VTKGeoMesh.h"
#include "ReadGmsh.h"
#include "CompMesh.h"
#include "Poisson.h"
#include "L2Projection.h"   // vai ser usado para função de contorno
#include "CompElement.h"
#include "GeoElement.h"
#include "Assemble.h"

using std::cout;
using std::endl;
using std::cin;

int main (){
    GeoMesh gmesh;          // ler a malha que criamos
        ReadGmsh read;
        read.Read(gmesh,"/home/emilia/Repositórios/FemCourseEigenClass2021/mainprograms/quads.msh");
        VTKGeoMesh plotmesh;
        plotmesh.PrintGMeshVTK(&gmesh, "quads.vtk");

    CompMesh cmesh(&gmesh);  // criar uma malha computacional
        MatrixDouble perm(3,3);
        perm.setZero();
        perm(0,0) = 1.;
        perm(1,1) = 1.;
        perm(2,2) = 1.;
        Poisson *mat1 = new Poisson(3, perm);

        MatrixDouble proj(1,1), val1(1,1), val2(1,1);
        proj.setZero();
        val1.setZero();
        val2.setOnes();
        L2Projection *bc_linha = new L2Projection(0, 2, proj, val1, val2);
        L2Projection *bc_point = new L2Projection(0, 1, proj, val1, val2);
        std::vector<MathStatement*>mathvec = {0, bc_point, bc_linha, mat1};
        cmesh.SetDefaultOrder(1);
        cmesh.SetMathVec(mathvec);
        cmesh.AutoBuild();
        cmesh.Resequence();

    //  Calculo da matriz de rigidez: mostrado em aula pelo professor
    for(auto cel:cmesh.GetElementVec())
    {   
        int nshape = cel->NShapeFunctions();
        auto nstate = cel->GetStatement()->NState();
        MatrixDouble ek(nstate * nshape, nstate * nshape);
        MatrixDouble ef(nstate * nshape, 1);
        ek.setZero();
        ef.setZero();

        auto gel = cel->GetGeoElement();
        auto nnodes = gel->NNodes();
        VecInt nodeindices;
        IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", "", "");
        IOFormat HeavyFmt(FullPrecision, 0, ", ", "\n", "{", "}", "{", "}");
        gel->GetNodes(nodeindices);
        std::cout << "element index " << cel->GetIndex() << std::endl;        
        std::cout << "coord = { ";
        for(auto in=0; in<nnodes; in++){
            GeoNode &node = gmesh.Node(nodeindices[in]);
            std::cout << "{ " << node.Co().format(CommaInitFmt) << "}";
            // std::cout << "{ " << node.Co() << "}";           
            if(in < nnodes-1) std::cout << ",";
        } 
        std::cout << "};\n";
        cel->CalcStiff(ek,ef);
        std::cout << "ek = " << ek.format(HeavyFmt) << ";\n";
        std::cout << "ef = " << ef.format(HeavyFmt) << ";\n";
        // std::cout << "ek = " << ek << ";\n";
        // std::cout << "ef = " << ef << ";";        
    }    

    plotmesh.PrintCMeshVTK(&cmesh, 2, "cmesh_quads.vtk"); 
    return 0;
}
