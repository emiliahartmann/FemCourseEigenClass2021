
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
    
    // std::cout << "order = " << order << " integral aproximada " << Integral <<
    // " erro " << 6.*M_PI-Integral << std::endl;


    GeoMesh gmesh;          // ler a malha que criamos
        ReadGmsh read;
        read.Read(gmesh,"/home/emilia/Repositórios/FemCourseEigenClass2021/mainprograms/malha.msh");
        VTKGeoMesh plotmesh;
        plotmesh.PrintGMeshVTK(&gmesh, "malha.vtk");

    CompMesh cmesh(&gmesh);  // criar uma malha computacional
        MatrixDouble perm(3,3);
        perm.setZero();
        perm(0,0) = 1.;
        perm(1,1) = 1.;
        perm(2,2) = 1.;
        Poisson *mat1 = new Poisson(1, perm);

        MatrixDouble proj(1,1), val1(1,1), val2(1,1);
        proj.setZero();
        val1.setZero();
        val2.setZero();
        // val2.setOnes();
        L2Projection *bc_linha = new L2Projection(0, 2, proj, val1, val2); // 0 -> 1 para condicao de Neumann
        L2Projection *bc_point = new L2Projection(0, 3, proj, val1, val2);
        std::vector<MathStatement*>mathvec (4);
        mathvec[0] = 0;
        mathvec[1] = mat1;
        mathvec[2] = bc_linha;
        mathvec[3] = bc_point;
        //= {0, mat1, bc_linha, bc_point}; // o meu material tem id igual a 1
        // cmesh.SetDefaultOrder(1); // se for 2, vai associar a funcoes quadraticas.
        cmesh.SetMathVec(mathvec); 
        cmesh.AutoBuild();
        cmesh.Resequence();
        cmesh.Solution() (0,0) = 1.;    // apenas para teste, excluir esta linha depois -- para visualizar as funcoes de forma
        plotmesh.PrintCMeshVTK(&cmesh, 2, "cmesh_malha.vtk"); 
        

        // Calculo da matriz de rigidez: mostrado em aula pelo professor
    // for(auto cel:cmesh.GetElementVec())
    // {
    //     MatrixDouble ek, ef;
    //     auto gel = cel->GetGeoElement();
    //     auto nnodes = gel->NNodes();
    //     VecInt nodeindices;
    //     // IOFormat CommaInitFmt(StreamPresicion, DontAlignCols, ", ", ", ", "", "", "", "");
    //     // IOFormat HeavyFmt(FullPresicion, 0, ", ", "\n", "{", "}", "{", "}");
    //     gel->GetNodes(nodeindices);
    //     std::cout << "element index " << cel->GetIndex() << std::endl;        
    //     std::cout << "coord = { ";
    //     for(auto in=0; in<nnodes; in++){
    //         GeoNode &node = gmesh.Node(nodeindices[in]);
    //         // std::cout << "{ " << node.Co().format(CommaInitFmt) << "}";
    //         std::cout << "{ " << node.Co() << "}";           
    //         if(in < nnodes-1) std::cout << ",";
    //     } 
    //     std::cout << "};\n";
    //     cel->CalcStiff(ek,ef);
    //     // std::cout << "ek = " << ek.format(HeavyFmt) << ";\n";
    //     // std::cout << "ef = " << ef.format(HeavyFmt) << ";\n";
    //     // std::cout << "ek = " << ek << ";\n";
    //     // std::cout << "ef = " << ef << ";";        
    // }    

    //  CalcStiff by Jefferson
    Analysis Analysis(&cmesh);
    Analysis.RunSimulation();


    // Assemblagem
    Assemble assemble(&cmesh);
    auto ne = assemble.NEquations();
    MatrixDouble globmat(ne, ne), rhs(ne, 1);
    assemble.Compute(globmat, rhs);
    return 0;

// +++++++++++++++++++++++++ To check the ek and ef values with mathematica +++++++++++++++++++++++++
// IntPointData data;
	
// 	data.axes.resize(2,3);
// 	data.axes.setZero();
// 	data.detjac = 1.;
// 	data.dphidksi.resize(2,3);
// 	data.dphidksi.setZero();
// 	data.dphidx.resize(2,3);
// 	data.gradx.resize(3,2);
// 	data.ksi.resize(2,1);
// 	data.phi.resize(3,1);
// 	data.weight = 1.;
// 	data.x.resize(3,1);
	
// 	// Inserir os valores correspondentes do Mathematica aqui: video 8 minutos
// 	data.phi[0]  = 0.3;
// 	data.phi[1]  = 0.3;
// 	data.phi[2]  = 0.4;
	
// 	data.dphidksi(0,0) = -1.;
// 	data.dphidksi(1,0) = -1.;		
// 	data.dphidksi(0,1) = 1.;
// 	data.dphidksi(1,1) = 0.;			
// 	data.dphidksi(0,2) = 0.;
// 	data.dphidksi(1,2) = 1.;	
	
// 	data.dphidx(0,0) = -1./M_SQRT2;
// 	data.dphidx(1,0) = 1./M_SQRT2-M_SQRT2;		
// 	data.dphidx(0,1) = 1./M_SQRT2;
// 	data.dphidx(1,1) = -1./M_SQRT2;		
// 	data.dphidx(0,2) = 0.;
// 	data.dphidx(1,2) = M_SQRT2;
			
// 	data.gradx(0,0) = 1.;		
// 	data.gradx(1,0) = 1.;
// 	data.gradx(0,1) = 0.;	
// 	data.gradx(1,1) = 1.;		
			
// 	data.axes(0,0) = 1./M_SQRT2;		
// 	data.axes(1,0) = 1./M_SQRT2;
// 	data.axes(0,1) = -1./M_SQRT2;	
// 	data.axes(1,1) = 1./M_SQRT2;	
	
// 	data.ksi[0] = 0.3;
// 	data.ksi[1] = 0.4;	
	
// 	data.weight = 0.2;
	
// 	data.x[0] = 0.3;
// 	data.x[1] = 0.7;			

// 	// Criar o objeto Poisson
// 	MatrixDouble perm(3,3);
// 	perm.setZero();
// 	perm(0,0) = 1.;
// 	perm(1,1) = 1.;
// 	perm(2,2) = 1.;
// 	Poisson matpoisson(1,perm);
// 	MatrixDouble ek(3,3), ef(3,1);
// 	ek.setZero();
// 	ef.setZero();
// 	matpoisson.Contribute(data, data.weight, ek, ef); // calculando a contribuicao da matriz de rigidez em um ponto somente da integracao
	
// 	std::cout << "ek/n" << ek << std::endl;
// 	std::cout << "ef/n" << ef << std::endl;
// 	return 0;
}
