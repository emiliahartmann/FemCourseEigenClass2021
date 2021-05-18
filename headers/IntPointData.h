//
//  IntPointData.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef IntPointData_h
#define IntPointData_h

#include "DataTypes.h"
#include "tpanic.h"

class IntPointData
{
public:
    
    VecDouble ksi;
    
    double weight;
    
    VecDouble phi;
    
    MatrixDouble dphidksi;
    
    VecDouble x;
    
    MatrixDouble gradx;
    
    MatrixDouble axes;
    
    double detjac;
    
    MatrixDouble dphidx;
    
    VecDouble solution;
    
    MatrixDouble dsoldksi;
    
    MatrixDouble dsoldx;
    
    VecDouble coefs;
    
    void ComputeSolution()
    {
        
        if(coefs.size()%phi.size())
        {
        //    DebugStop();
        }
        int nstate = coefs.size()/phi.size();
        solution.resize(nstate);
        solution.setZero();
        dsoldksi.resize(dphidx.rows(), nstate);
        dsoldx.resize(dphidx.rows(), nstate);
        dsoldx.setZero();
        dsoldksi.setZero();
         
        int dim = dphidx.rows();
        for (int iphi=0; iphi<phi.size(); iphi++) {
            double phival = phi[iphi];
            for (int istate=0; istate<nstate; istate++) {
                solution[istate] += phival*coefs[iphi*nstate+istate];
                for (int d=0; d < dim; d++) {
                    dsoldksi(d,istate) += coefs[iphi*nstate+istate]*dphidksi(d,iphi);
                    dsoldx(d,istate) += coefs[iphi*nstate+istate]*dphidx(d,iphi);
                }
            }
        }
    }
    
};
#endif /* IntPointData_h */
