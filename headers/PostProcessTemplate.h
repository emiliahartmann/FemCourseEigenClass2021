//
//  PostProcessTemplate.h
//  FemCourse
//
//  Created by Philippe Devloo on 15/05/18.
//

#ifndef PostProcessTemplate_h
#define PostProcessTemplate_h

#include "MathStatement.h"
#include "Poisson.h"
#include "Analysis.h"
#include "L2Projection.h"
#include "IntPointData.h"
#include "PostProcess.h"
#include <list>

class PostProcess;

template<class math>
class PostProcessTemplate: public PostProcess
{

    protected:
    
    std::vector<typename math::PostProcVar> scalarvariables;
    std::vector<typename math::PostProcVar> vectorvariables;
    
    std::vector<std::string> vectornames;
    std::vector<std::string> scalarnames;

    
    public:
    
    PostProcessTemplate() : PostProcess(){
    }
    
    PostProcessTemplate(const PostProcessTemplate &copy) : PostProcess(copy) {
    }
    
    ~PostProcessTemplate(){
        
    }
    
    
    PostProcessTemplate &operator=(const PostProcessTemplate &cp){
        return *this;
    }
    
    PostProcessTemplate(Analysis *Ref) : PostProcess(Ref){
    }
    
    
    virtual std::vector<std::string> Vectornames() const{
        return vectornames;
    }

    virtual std::vector<std::string> Scalarnames() const{
        return scalarnames;
    }
    
    virtual VecInt VectorvariablesIds() const{
        math Statement;
        VecInt VecVar(NumVectorVariables(),0);
        
        for (int index = 0; index< NumVectorVariables(); index++) {
            VecVar[index]=Statement.VariableIndex(vectorvariables[index]);
        }
        
        return VecVar;
    }

    virtual VecInt ScalarvariablesIds() const{
        math Statement;
        VecInt ScalVar(NumScalarVariables(),0);
        
        for (int index = 0; index< NumScalarVariables(); index++) {
            ScalVar[index]=Statement.VariableIndex(scalarvariables[index]);
        }
        
        return ScalVar;
    }
    
    virtual void AppendVariable(const std::string &name){
        math Statement;
        typename math::PostProcVar var = Statement.VariableIndex(name);
        AppendVariable(var);
    
        int nsol = Statement.NSolutionVariables(var);
        
        if (nsol==1) {
            scalarnames.push_back(name);
        }else{
            vectornames.push_back(name);
        }
        
    }
    
    
    void AppendVariable(typename math::PostProcVar var){
        math Statement;
        int nsol = Statement.NSolutionVariables(var);
        
        if (nsol==1) {
            scalarvariables.push_back(var);
        }else{
            vectorvariables.push_back(var);
        }
        
    }
    
    VecDouble PostProcResult(MathStatement &mathStatement, unsigned int varIndex, const IntPointData &data) const {
        math *locptr = dynamic_cast<math *> (&mathStatement);
        
        if(!locptr) DebugStop();
        const int numScalarVariables = NumScalarVariables();
        VecDouble sol;
        locptr->PostProcessSolution(data, varIndex < numScalarVariables? scalarvariables[varIndex] : vectorvariables[varIndex-numScalarVariables],sol);
        return sol;
    }
    
    virtual inline unsigned int NumScalarVariables() const {
        return scalarvariables.size();
    }
    
    virtual  inline unsigned int NumVectorVariables() const {
        return vectorvariables.size();
    }
    
    virtual inline unsigned int NumVariables() const {
        return NumScalarVariables() + NumVectorVariables();
    }
};
#endif /* PostProcessTemplate_h */

template class PostProcessTemplate<Poisson>;
template class PostProcessTemplate<L2Projection>;
