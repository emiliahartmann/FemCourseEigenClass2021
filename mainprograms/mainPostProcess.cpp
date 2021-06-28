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
#include "PostProcessTemplate.h"
#include "PostProcess.h"

using std::cout;
using std::endl;
using std::cin;

int main()
{
	IntPointData data;
	
	data.axes.resize(2,3);
	data.axes.setZero();
	data.detjac = 1.;
	data.dphidksi.resize(2,3);
	data.dphidksi.setZero();
	data.dphidx.resize(2,3);
	data.gradx.resize(3,2);
	data.ksi.resize(2,1);
	data.phi.resize(3,1);
	data.weight = 1.;
	data.x.resize(3,1);
	
	// Inserir os valores correspondentes do Mathematica aqui: video 8 minutos
	data.phi[0]  = 0.3;
	data.phi[1]  = 0.3;
	data.phi[2]  = 0.4;
	
	data.dphidksi(0,0) = -1.;
	data.dphidksi(1,0) = -1.;		
	data.dphidksi(0,1) = 1.;
	data.dphidksi(1,1) = 0.;			
	data.dphidksi(0,2) = 0.;
	data.dphidksi(1,2) = 1.;	
	
	data.dphidx(0,0) = -1./M_SQRT2;
	data.dphidx(1,0) = 1./M_SQRT2-M_SQRT2;		
	data.dphidx(0,1) = 1./M_SQRT2;
	data.dphidx(1,1) = -1./M_SQRT2;		
	data.dphidx(0,2) = 0.;
	data.dphidx(1,2) = M_SQRT2;
			
	data.gradx(0,0) = 1.;		
	data.gradx(1,0) = 1.;
	data.gradx(0,1) = 0.;	
	data.gradx(1,1) = 1.;		
			
	data.axes(0,0) = 1./M_SQRT2;		
	data.axes(1,0) = 1./M_SQRT2;
	data.axes(0,1) = -1./M_SQRT2;	
	data.axes(1,1) = 1./M_SQRT2;	
	
	data.ksi[0] = 0.3;
	data.ksi[1] = 0.4;	
	
	data.weight = 0.2;
	
	data.x[0] = 0.3;
	data.x[1] = 0.7;			

data.solution.resize(1);
data.solution[0] = 1;
data.dsoldx.resize(2,1);
data.dsoldx(0,0) = -5.;
data.dsoldx(1,0) = 2.;

	// Criar o objeto Poisson
	MatrixDouble perm(3,3);
	perm.setZero();
	perm(0,0) = 1.;
	perm(1,1) = 1.;
	perm(2,2) = 1.;
	Poisson matpoisson(1,perm);

	PostProcessTemplate<Poisson> postprocess;   // postprocess e o objeto

    VecDouble Solout;	// era solout
    postprocess.AppendVariable("Flux");
    postprocess.AppendVariable("Sol");
    postprocess.AppendVariable("DSol");
    postprocess.AppendVariable("SolExact");
    postprocess.AppendVariable("Force");
    postprocess.AppendVariable("DSolExact");

    int64_t nscal = postprocess.NumScalarVariables();
    int64_t nvecs = postprocess.NumVectorVariables();

	std::cout << "Scalar post processing/n";
    for (int i = 0; i < nscal; i++){
        auto name = postprocess.Scalarnames()[i];
        std::cout << "name " << name << " post processed " << postprocess.PostProcResult(matpoisson, i , data) << std::endl;
    }

    	std::cout << "Vector post processing/n";
    for (int i = 0; i < nvecs; i++){
        auto name = postprocess.Vectornames()[i];
        std::cout << "name " << name << " post processed " << postprocess.PostProcResult(matpoisson, i+nscal , data) << std::endl;
    }

	return 0;
}