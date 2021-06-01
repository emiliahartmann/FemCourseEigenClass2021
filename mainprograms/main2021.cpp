

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

using std::cout;
using std::endl;
using std::cin;

int main (){

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

    return 0;
}
