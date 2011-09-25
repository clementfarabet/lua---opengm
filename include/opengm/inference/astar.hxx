/// OpenGM. Copyright (c) 2010 by Bjoern Andres and Joerg Hendrik Kappes.
///
/// This software was developed by Bjoern Andres and Joerg Hendrik Kappes.
/// Enquiries shall be directed to:
/// bjoern.andres@iwr.uni-heidelberg.de, kappes@math.uni-heidelberg.de
///
/// Author(s) of this file: Joerg Hendrik Kappes
///
/// All advertising materials mentioning features or use of this software must
/// display the following acknowledgement: ``This product includes the OpenGM
/// library developed by Bjoern Andres and Joerg Hendrik Kappes. Please direct 
/// enquiries concerning OpenGM to bjoern.andres@iwr.uni-heidelberg.de,
/// kappes@math.uni-heidelberg.de''.
///
/// Redistribution and use in source and binary forms, with or without
/// modification, are permitted provided that the following conditions are met:
///
/// - Redistributions of source code must retain the above copyright notice,
///   this list of conditions and the following disclaimer.
/// - Redistributions in binary form must reproduce the above copyright notice, 
///   this list of conditions and the following disclaimer in the documentation
///   and/or other materials provided with the distribution.
/// - All advertising materials mentioning features or use of this software must 
///   display the following acknowledgement: ``This product includes the OpenGM
///   library developed by Bjoern Andres and Joerg Hendrik Kappes. Please direct 
///   enquiries concerning OpenGM to bjoern.andres@iwr.uni-heidelberg.de,
///   kappes@math.uni-heidelberg.de''.
/// - The names of the authors must not be used to endorse or promote products 
///   derived from this software without specific prior written permission.
///
/// THIS SOFTWARE IS PROVIDED BY THE AUTHORS ``AS IS'' AND ANY EXPRESS OR IMPLIED 
/// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
/// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
/// EVENT SHALL THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
/// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
/// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
/// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
/// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
/// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/// 
#pragma once
#ifndef OPENGM_ASTAR_HXX
#define OPENGM_ASTAR_HXX

#include <cmath>
#include <vector>
#include <list>
#include <set>
#include <iostream>
#include <algorithm> 
#include <iostream> 
#include <functional>

#include "opengm/opengm.hxx"
#include "opengm/maxdistance.hxx"
#include "opengm/graphicalmodel.hxx"
#include "opengm/exitcondition.hxx"
#include "opengm/inference/inference.hxx"
#include "opengm/inference/beliefpropagation.hxx"
#include <opengm/timing.hxx>

namespace opengm {

// node of the search tree for the a-star search
template<class Factor> struct AStarNode {
    typename std::vector<size_t>    conf;
    typename Factor::value_type     value;
};

// a-star search
template<class GM,class ACC>
class AStar : public Inference<GM,ACC>
{
public:
    typedef GM                                  gm_type;
    typedef ACC                                 Accumulation;
    typedef typename gm_type::factor_type       factor_type;
    typedef typename gm_type::factor_type       Factor;
    typedef typename factor_type::value_type    value_type;   
    typedef typename factor_type::space_type    space_type;
    typedef typename space_type::state_type     state_type;
    typedef typename gm_type::Operator          Operator;
    typedef typename std::vector<state_type>    ConfVec ;  
    typedef typename ConfVec::iterator          ConfVecIt;

    struct Parameter {
        Parameter()
            {
                maxHeapSize_    = 3000000;
                maxTimeMs_      = std::numeric_limits<double>::infinity();
                numberOfOpt_    = 1;
                objectiveBound_ = Accumulation::template neutral<value_type>();
                heuristic_      = Parameter::DEFAULTHEURISTIC; 
            };
        void addTreeFactorId(size_t id)
            { treeFactorIds_.push_back(id); }
 
        static const size_t DEFAULTHEURISTIC = 0;
        static const size_t FASTHEURISTIC = 1;
        static const size_t STANDARDHEURISTIC = 2;
        size_t maxHeapSize_;
        std::vector<size_t> nodeOrder_;
        std::vector<size_t> treeFactorIds_;
        double              maxTimeMs_;
        size_t              numberOfOpt_;
        value_type          objectiveBound_;
        //ExitConditions<value_type,Operator,Accumulation> exitCond_;
        size_t heuristic_;
    };

    AStar(const GM& gm, Parameter para = Parameter());
    virtual std::string name() const {return "AStar";} 
    const gm_type& graphicalModel() const;

    virtual InferenceTermination infer();
    template<class VisitorType> InferenceTermination infer(VisitorType& vistitor);
    virtual InferenceTermination marginal(const size_t&,Factor& out)const        {return UNKNOWN;}
    virtual InferenceTermination factorMarginal(const size_t&, Factor& out)const {return UNKNOWN;}
    virtual InferenceTermination arg(std::vector<state_type>& v, const size_t& = 1)const;
    virtual InferenceTermination args(std::vector< std::vector<state_type> >& v)const;

private:
    const GM&                         gm_;
    Parameter                         parameter_;
    std::vector<AStarNode<Factor> >   array_;

    std::vector<size_t>               numStates_;
    size_t                            numNodes_;
    std::vector<Factor>               treeFactor_;
    std::vector<Factor>               optimizedFactor_;
    std::vector<ConfVec >             optConf_;
    std::vector<bool>                 isTreeFactor_;
    value_type                        bound2_;

    template<class VisitorType> void  expand(VisitorType& vistitor);
    std::vector<value_type>           fastHeuristic(ConfVec conf);

    inline static bool                comp1(const AStarNode<Factor>& a, const AStarNode<Factor>& b) 
        {return  Accumulation::ibop(a.value,b.value);};
    inline static bool                comp2(const AStarNode<Factor>& a, const AStarNode<Factor>& b)  
        { return  Accumulation::bop(a.value,b.value);}; 
    //inline static value_type          operation(value_type a, value_type b){value_type t; Operator::op(a,b,t);return t;};
    inline static value_type          better(value_type a, value_type b)   {return Accumulation::op(a,b);};
    inline static value_type          wrose(value_type a,  value_type b)   {return Accumulation::iop(a,b);};
};

// visitors
/*
template<class AStar>
class AStarVisitor {
public:
  typedef AStar astar_type;
  typedef typename astar_type::value_type value_type;
 
  AStarVisitor();
  void operator()(const astar_type&, const std::vector<size_t>& conf, const size_t& heapsize, const value_type&, const double& runtime) const;
};
*/
  template<class AStar, bool Verbose=false>
class AStarVisitor {
public:
  typedef AStar astar_type;
  typedef typename astar_type::value_type value_type;

  AStarVisitor();
  void operator()(const astar_type&, const std::vector<size_t>& conf, const size_t& heapsize, const value_type& bound1, const value_type& bound2, const double& runtime);
private:
  size_t step_;
};



//*****************
//** DEFINITIONS **
//*****************

template<class GM, class ACC >
AStar<GM,ACC>
    ::AStar
    (
    const GM& gm,
    Parameter para
    ):gm_(gm)
{
    parameter_ = para;
    if( parameter_.heuristic_ == Parameter::DEFAULTHEURISTIC){
      if(gm_.factorOrder()<=2)
	parameter_.heuristic_ = Parameter::FASTHEURISTIC;
      else
	parameter_.heuristic_ = Parameter::STANDARDHEURISTIC;
    }
    assert(parameter_.heuristic_ == Parameter::FASTHEURISTIC ||	parameter_.heuristic_ == Parameter::STANDARDHEURISTIC );

    bound2_ = std::numeric_limits<value_type>::infinity();

    //Set variables
    isTreeFactor_.resize(gm_.numberOfFactors());
    numStates_.resize(gm_.space().dimension());
    numNodes_ = gm_.space().dimension();
    for(size_t i=0; i<numNodes_;++i)
        numStates_[i] = gm_.space().numberOfStates(i);

    //Check nodeOrder
    if(parameter_.nodeOrder_.size()==0){
        parameter_.nodeOrder_.resize(numNodes_);
        for(size_t i=0; i<numNodes_; ++i)
            parameter_.nodeOrder_[i]=i;
    }
    if(parameter_.nodeOrder_.size()!=numNodes_)
        throw std::runtime_error("Node order does not fit to the model");

    assert(std::set<size_t>(parameter_.nodeOrder_.begin(), parameter_.nodeOrder_.end()).size()==numNodes_);
    for(size_t i=0;i<numNodes_; ++i){
        assert(parameter_.nodeOrder_[i]<numNodes_);
        assert(parameter_.nodeOrder_[i]>=0);
    }

    //Check FactorIds
    if(parameter_.treeFactorIds_.size()==0){
      //Select tree factors
      for(size_t i=0; i<gm_.numberOfFactors(); ++i){
	if((gm_[i].numberOfVariables()==2) && 
	   (gm_[i].variableIndex(0)==parameter_.nodeOrder_.back() || gm_[i].variableIndex(1)==parameter_.nodeOrder_.back())
	  )
	   parameter_.addTreeFactorId(i);
      }
    }
    for(size_t i=0; i<parameter_.treeFactorIds_.size(); ++i)
        assert(gm_.numberOfFactors() > parameter_.treeFactorIds_[i]);

    //compute optimed factors
    optimizedFactor_.resize(gm_.numberOfFactors());
    for(size_t i=0; i<gm_.numberOfFactors(); ++i){
        if(gm_[i].numberOfVariables()<=1) continue; 
        std::vector<size_t> index(gm_[i].numberOfVariables()); 
        gm_[i].variableIndices(index.begin());
        optimizedFactor_[i].assign(gm_.space(),index.end()-1, index.end());
        gm_[i].template accumulate<ACC>(index.begin()+1,index.end(),optimizedFactor_[i]);
        assert(optimizedFactor_[i].numberOfVariables()==1);
        assert(optimizedFactor_[i].variableIndex(0)==index[0]);

    }
    //PUSH EMPTY CONFIGURATION TO HEAP
    AStarNode<Factor> a;
    a.conf.resize(0); 
    a.value = 0;      
    array_.push_back(a); 
    make_heap(array_.begin(), array_.end(), comp1);

    //Check if maximal order is smaller equal 2, otherwise fall back to naive computation of heuristic
    if(parameter_.heuristic_ == parameter_.FASTHEURISTIC){
        for(size_t i=0; i<parameter_.treeFactorIds_.size(); ++i){
            if(gm_[parameter_.treeFactorIds_[i]].numberOfVariables()>2)
                throw std::runtime_error("Heuristic includes factor with order>2");
        }
    }

    //Init treefactor structure
    treeFactor_.clear();
    for(size_t i=0; i<gm_.numberOfFactors(); ++i)
        isTreeFactor_[i] = false;
    for(size_t i=0; i<parameter_.treeFactorIds_.size(); ++i){
        int factorId = parameter_.treeFactorIds_[i];
        isTreeFactor_[factorId] = true;
        treeFactor_.push_back(gm_[factorId]);
    }
}

///////////////////////////////////////
///////////////////////////////////////

template <class GM, class ACC>
InferenceTermination 
AStar<GM,ACC>::infer()
{
    AStarVisitor<AStar<gm_type, ACC> > v;
    return infer(v);
}

template<class GM, class ACC>
template<class VisitorType>
InferenceTermination AStar<GM,ACC>::infer(VisitorType& visitor)
{ 
  double time=0;
  USETICTOC;
  TIC;
  optConf_.resize(0);
  while(array_.size()>0){
    if(parameter_.numberOfOpt_ == optConf_.size()){
      return NORMAL;
    }
    while(array_.front().conf.size() < numNodes_){
      expand(visitor);
     
      time += TOCN;     
      visitor(*this, array_.front().conf, array_.size(), array_.front().value, bound2_, time);
      TIC;
      if(time>parameter_.maxTimeMs_)
      	return TIMEOUT; 
    }
    
    value_type  value = array_.front().value;
       
    std::vector<size_t> conf(numNodes_);
    for(size_t n=0; n<numNodes_; ++n){
      conf[parameter_.nodeOrder_[n]] = array_.front().conf[n];
    }

    if(ACC::bop(parameter_.objectiveBound_, value)){
      return NORMAL;
    }
    optConf_.push_back(conf); 
    pop_heap(array_.begin(), array_.end(),  comp1); //greater<Factor,Accumulation>);
    array_.pop_back();
  }
  return UNKNOWN;
}

template<class GM, class ACC>
InferenceTermination AStar<GM, ACC>
    ::arg(ConfVec& conf, const size_t& n)const
{
        if(n>optConf_.size()){
            conf.resize(0);
            return UNKNOWN;
        }
        //conf.resize(opt_conf[n-1].size());
        conf=optConf_[n-1];
        return NORMAL;
}

template<class GM, class ACC>
InferenceTermination AStar<GM,ACC>
    ::args(std::vector<std::vector<state_type> >& conf)const
{
        conf=optConf_;
        return NORMAL;
}

template<class GM, class ACC>
template<class VisitorType>
void AStar<GM, ACC>::expand(VisitorType& visitor)
{ 
   //CHECK HEAP SIZE
    if(array_.size()>parameter_.maxHeapSize_*0.99){ 
        partial_sort(array_.begin(), array_.begin()+(int)(parameter_.maxHeapSize_/2), array_.end(),  comp2);
        array_.resize((int)(parameter_.maxHeapSize_/2));
    }
   
    //GET HEAP HEAD
    AStarNode<Factor> a           = array_.front();
    size_t            subconfsize = a.conf.size();
 
    //REMOVE HEAD FROM HEAP
    assert(array_.size()>0);
    pop_heap(array_.begin(), array_.end(),  comp1); //greater<Factor,Accumulation>);
    array_.pop_back();
     
    if( parameter_.heuristic_ == parameter_.STANDARDHEURISTIC){
        //BUILD GRAPHICAL MODEL FOR HEURISTC CALCULATION
        GM tgm = gm_;
        std::vector<size_t> variableIndices(subconfsize);
        std::vector<size_t> values(subconfsize);

        for(size_t i =0; i<subconfsize ; ++i){
            variableIndices[i] = parameter_.nodeOrder_[i];
            values[i]          = a.conf[i];
        }
        tgm.introduceEvidence(variableIndices.begin(), variableIndices.end(),values.begin());    
        std::vector<size_t> varInd;
        for(size_t i=0; i<tgm.numberOfFactors(); ++i){

            //treefactors will not be modyfied
            if(isTreeFactor_[i]) continue;

            varInd.clear();
            size_t nvar = tgm[i].numberOfVariables(); 

            //factors depending from 1 variable can be include to the tree
            if(nvar<=1) continue;

	    tgm[i] = optimizedFactor_[i];
        } 

        //double c = 0; //FIXME unused
        BeliefPropagation<GM, ACC, opengm::MaxDistance> bp(tgm);
        bp.inferAsynchronous();

	if(true){
	  std::vector<size_t> conf(numNodes_);
	  bp.arg(conf);
	  for(size_t i =0; i<subconfsize ; ++i){
            conf[i]          = a.conf[i];
	  }
	  bound2_= ACC::op(gm_.evaluate(conf),bound2_);
	}


        std::vector<size_t> conf(numNodes_);

        a.conf.resize(subconfsize+1); 
        for(size_t i=0; i<numStates_[subconfsize]; ++i){
            a.conf[subconfsize] = i;
            bp.constrainedOptimum(parameter_.nodeOrder_,a.conf,conf);
            a.value   = tgm.evaluate(conf);
            array_.push_back(a);    
            push_heap(array_.begin(), array_.end(),  comp1); //greater<Factor,Accumulation>) ;
        } 
    }

    if( parameter_.heuristic_ == parameter_.FASTHEURISTIC){
        std::vector<size_t> conf(subconfsize);
        for(size_t i=0;i<subconfsize;++i)
            conf[i] = a.conf[i];   
        std::vector<value_type> bound = fastHeuristic(conf);
        a.conf.resize(subconfsize+1);
        for(size_t i=0; i<numStates_[parameter_.nodeOrder_[subconfsize]]; ++i){
            a.conf[subconfsize] = i;
            a.value             = bound[i];
            //if(bound[i]<10){
            array_.push_back(a);   
            push_heap(array_.begin(), array_.end(),  comp1); //greater<Factor,Accumulation>) ;
            //}
        }

    }


}

template<class GM, class ACC>  
std::vector<typename AStar<GM, ACC>::value_type> 
AStar<GM, ACC>::fastHeuristic(typename AStar<GM, ACC>::ConfVec conf)
{ 
    std::list<size_t>                 factorList;
    std::vector<size_t>               nodeDegree(numNodes_,0);
    std::vector<int>                  nodeLabel(numNodes_,-1);
    std::vector<std::vector<value_type > > nodeEnergy(numNodes_); 
    size_t                            nextNode = parameter_.nodeOrder_[conf.size()];

    for(size_t i=0; i<numNodes_; ++i){
        nodeEnergy[i].resize(numStates_[i]); //the energy passed to node i
        for(size_t j=0;j<numStates_[i];++j)
            Operator::neutral(nodeEnergy[i][j]); 
    }    
    for(size_t i=0;i<conf.size();++i){
        nodeLabel[parameter_.nodeOrder_[i]] = conf[i];
    }

    //First run:
    // * add unarry function
    // * add pairwise functions with at least one observed node
    // * add the approximation for pairwise none-tree edges
    for(size_t i=0; i<gm_.numberOfFactors(); ++i){

        Factor f    = gm_[i];
        size_t nvar = f.numberOfVariables(); 

        //factors depending from 1 variable can be include 
        if(nvar==1){
            int index = f.variableIndex(0);
            if(nodeLabel[index]>=0){
                nodeEnergy[index].resize(1);
                //nodeEnergy[index][0] = operatipon(f(nodeLabel[index]), nodeEnergy[index][0]);
                Operator::op(f(nodeLabel[index]), nodeEnergy[index][0]);
            }
            else{  
                assert(numStates_[index]==nodeEnergy[index].size());
                for(size_t j=0;j<numStates_[index];++j){
                    //nodeEnergy[index][j] = operation(f(j),nodeEnergy[index][j]);
                    Operator::op(f(j),nodeEnergy[index][j]);
                }
            } 
        }
        if(nvar==2){
            size_t index1 = f.variableIndex(0);
            size_t index2 = f.variableIndex(1);

            if(nodeLabel[index1]>=0){
                if(nodeLabel[index2]>=0){
                    nodeEnergy[index1].resize(1);
                    //nodeEnergy[index1][0] = operation(f(nodeLabel[index1],nodeLabel[index2]),nodeEnergy[index1][0]); 
                    Operator::op(f(nodeLabel[index1],nodeLabel[index2]),nodeEnergy[index1][0]); 
                }
                else{
                    assert(numStates_[index2]== nodeEnergy[index2].size());
                    for(size_t j=0;j<numStates_[index2];++j){
                        //nodeEnergy[index2][j] = operation(f(nodeLabel[index1],j), nodeEnergy[index2][j]);
                        Operator::op(f(nodeLabel[index1],j), nodeEnergy[index2][j]);
                    }
                }
            }
            else if(nodeLabel[index2]>=0){
                assert(numStates_[index1]==nodeEnergy[index1].size());
                for(size_t j=0;j<numStates_[index1];++j){
                    //nodeEnergy[index1][j] = operation(f(j,nodeLabel[index2]),nodeEnergy[index1][j]);
                    Operator::op(f(j,nodeLabel[index2]),nodeEnergy[index1][j]);
                }
            }	
            else if(isTreeFactor_[i]){
                factorList.push_front(i);
                ++nodeDegree[index1];
                ++nodeDegree[index2];
                continue;
            }
            else{	    
                for(size_t j=0;j<numStates_[index1];++j){
                    //nodeEnergy[index1][j] = operation(optimizedFactor_[i](j), nodeEnergy[index1][j]);
                    Operator::op(optimizedFactor_[i](j), nodeEnergy[index1][j]);
                }
            }
        }
        if(nvar>2){
            bool covered = true;
            std::vector<size_t> state(nvar);
            for(size_t j=0; j<nvar; ++j){
                if(nodeLabel[f.variableIndex(j)]<0){
                    state[j] = nodeLabel[f.variableIndex(j)];
                    covered = false;
                }
            }
            if(covered)
                nodeEnergy[f.variableIndex(0)][0] = f(state.begin());
            else{
                for(size_t j=0;j<numStates_[f.variableIndex(0)];++j){
                    //nodeEnergy[f.variableIndex(0)][j] = operation(optimizedFactor_[i](j), nodeEnergy[f.variableIndex(0)][j]);
                    Operator::op(optimizedFactor_[i](j), nodeEnergy[f.variableIndex(0)][j]);
                }
            }
        }
    } 

    nodeDegree[nextNode] += numNodes_;

    // Start dynamic programming to solve the treestructured problem.
    while(factorList.size()>0){
        size_t    id  = factorList.front();
        factorList.pop_front();
        Factor f      = gm_[id];
        size_t    index1 = f.variableIndex(0);
        size_t    index2 = f.variableIndex(1);
        typename Factor::value_type temp;

        assert(index1<numNodes_);
        assert(index2<numNodes_);
        assert(f.space().numberOfStates(index1) == numStates_[index1]);
        assert(f.space().numberOfStates(index1) == numStates_[index1]);


        if(nodeDegree[index1]==1){
            typename Factor::value_type min;
            assert(numStates_[index2]==nodeEnergy[index2].size());
            for(size_t j2=0;j2<numStates_[index2];++j2){
                ACC::neutral(min);
                assert(numStates_[index1] == nodeEnergy[index1].size());
                for(size_t j1=0;j1<numStates_[index1];++j1){
                    Operator::op(f(j1,j2),nodeEnergy[index1][j1],temp);
                    ACC::op(min,temp,min);
                }
                //nodeEnergy[index2][j2] = operation(min,nodeEnergy[index2][j2]);
                Operator::op(min,nodeEnergy[index2][j2]); 
            } 
            --nodeDegree[index1];
            --nodeDegree[index2];
            nodeEnergy[index1].resize(1);
            Operator::neutral(nodeEnergy[index1][0]);
        }
        else if(nodeDegree[index2]==1){
            typename Factor::value_type min;
            assert(numStates_[index1]==nodeEnergy[index1].size());
            for(size_t j1=0;j1<numStates_[index1];++j1){
                ACC::neutral(min);
                assert(numStates_[index2]==nodeEnergy[index2].size());
                for(size_t j2=0;j2<numStates_[index2];++j2){
                    Operator::op(f(j1,j2),nodeEnergy[index2][j2],temp);
                    ACC::op(min,temp,min);
                    //if(min>f(j1,j2)*node_energy[index2][j2]) min=f(j1,j2)*node_energy[index2][j2];
                }
                //nodeEnergy[index1][j1] = operation(min,nodeEnergy[index1][j1]);
                Operator::op(min,nodeEnergy[index1][j1]);  
            }
            --nodeDegree[index1];
            --nodeDegree[index2];
            nodeEnergy[index2].resize(1);
            Operator::neutral(nodeEnergy[index2][0]);
        }
        else{
            factorList.push_back(id);
        }
    }

    //Evaluate
    value_type result;     
    value_type min;
    Operator::neutral(result);
    std::vector<value_type > bound;
    for(size_t i=0;i<numNodes_;++i){
        if(i==nextNode) continue;
        ACC::neutral(min);
        for(size_t j=0; j<nodeEnergy[i].size();++j)
            ACC::op(min,nodeEnergy[i][j],min);
        //result = operation(result,min);
        Operator::op(min,result);
    }    

    bound.resize(nodeEnergy[nextNode].size());
    for(size_t j=0; j<nodeEnergy[nextNode].size();++j){
        //bound[j] = operation(nodeEnergy[nextNode][j],result);	
        Operator::op(nodeEnergy[nextNode][j],result,bound[j]);
    }    
    return bound;
} 

template<class GM, class ACC>
inline const typename AStar<GM, ACC>::gm_type&
AStar<GM, ACC>::graphicalModel() const
{
    return gm_;
}  


// implementation of visitors
/*
template<class AStar>
inline AStarVisitor<AStar>::AStarVisitor()
{ }

template<class AStar>
inline void
AStarVisitor<AStar>::operator()
(
    const typename AStarVisitor<AStar>::astar_type& astar,
    const std::vector<size_t>& conf,
    const size_t& heapsize,
    const typename AStarVisitor<AStar>::value_type& bound,
    const double& runtime
) const
{ }

*/

template<class AStar, bool Verbose>
inline AStarVisitor<AStar,Verbose>::AStarVisitor()
: step_(0)
{ 
}


template<class AStar, bool Verbose>
inline void
AStarVisitor<AStar,Verbose>::operator() 
(
 const typename AStarVisitor<AStar,Verbose>::astar_type& astar,
    const std::vector<size_t>& conf,
    const size_t& heapsize,
    const typename AStarVisitor<AStar,Verbose>::value_type& bound1,
    const typename AStarVisitor<AStar,Verbose>::value_type& bound2,
    const double& runtime
) 
{
  if(Verbose){
    ++step_;
    std::cout << "step :" << step_
              << "    time = " << runtime/1000.0 << "s"  
	      << "    heapsize = "<< heapsize 
	      << "    " << bound1 << " <= E(x) <= " << bound2 
	      << std::endl;
  }
}

} // namespace opengm

#endif // #ifndef OPENGM_ASTAR_HXX

