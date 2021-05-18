//
//  Analysis.h
//  FemCourse
//
//  Created by Philippe Devloo on 08/05/18.
//

#ifndef Analysis_h
#define Analysis_h


#include "DataTypes.h"
#include <string>
#include <iostream>

//#include <armadillo>
//using namespace arma;

class CompMesh;

class PostProcess;

class Analysis
{
protected:
    CompMesh *cmesh;
    
    MatrixDouble Solution;
    
    MatrixDouble GlobalSystem;
    
    MatrixDouble RightHandSide;
    
public:
    
    Analysis();
    
    Analysis(const Analysis &cp);
    
    Analysis &operator=(const Analysis &cp);
    
    ~Analysis();
    
    Analysis(CompMesh *cmesh);
    
    void SetMesh(CompMesh *cmesh);
    
    CompMesh *Mesh() const;
    
    void RunSimulation();
    
    void PostProcessSolution(const std::string &filename, PostProcess &defPostProc) const;
    
    VecDouble PostProcessError(std::ostream &out, PostProcess &defPostProc) const;
    
};

#endif /* Analysis_h */
