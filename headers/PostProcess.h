//
//  PostProcess.h
//  FemCourse
//
//  Created by Philippe Devloo on 15/05/18.
//

#ifndef PostProcess_h
#define PostProcess_h
#include <functional>
#include "DataTypes.h"

class Analysis;

class PostProcess
{
protected:
    
    Analysis *Reference;

    // Pointer to Exact solution function, it is necessary to calculating errors
    std::function<void (const VecDouble &loc, VecDouble &result, MatrixDouble &deriv)> fExact;
    
public:
    
    PostProcess(){
        Reference=0;
    }

    PostProcess(const PostProcess &copy){
        Reference=copy.Reference;
    }
    
    ~PostProcess(){
        
    }
    
    
    PostProcess &operator=(const PostProcess &cp){
        Reference=cp.Reference;
        return *this;
    }
    
    PostProcess(Analysis *Ref){
        Reference=Ref;
    }
    
    virtual void Write(std::string filename){
        
    }
    
    virtual void SetExact(std::function<void (const VecDouble &loc, VecDouble &result, MatrixDouble &deriv)> Exact){
        fExact=Exact;
    }
    
    virtual std::function<void (const VecDouble &loc, VecDouble &result, MatrixDouble &deriv)> GetExact(){
        return fExact;
    }
    
    virtual std::vector<std::string> Vectornames() const = 0;
    
    virtual std::vector<std::string> Scalarnames() const = 0;
    
    virtual VecInt VectorvariablesIds() const = 0;
    
    virtual VecInt ScalarvariablesIds() const = 0;
    
    virtual void AppendVariable(const std::string &name) = 0;
    
    virtual inline unsigned int NumScalarVariables() const = 0;
    
    virtual  inline unsigned int NumVectorVariables() const = 0;
    
    virtual inline unsigned int NumVariables() const = 0;
    
    
};
#endif /* PostProcess_h */
