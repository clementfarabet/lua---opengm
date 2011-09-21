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
#ifndef OPENGM_TREEREWEIGTHEDBELIEFPROPAGATION_HXX
#define OPENGM_TREEREWEIGTHEDBELIEFPROPAGATION_HXX

#include <vector>
#include <map>
#include <list>
#include <set>

#include "opengm/opengm.hxx"
#include "opengm/graphicalmodel.hxx"
#include "opengm/factorbuffer.hxx"
#include "opengm/inference/inference.hxx"
#include "opengm/decomposer.hxx"

namespace opengm {

/// Wrapper class for variable nodes
template<class Factor, class OP, class ACC>
class VariableHullTRBP {
public:
    typedef typename Factor::space_type space_type;
    typedef typename Factor::value_type value_type;

    VariableHullTRBP();
    void assign(const space_type&, const size_t&, const std::vector<size_t>&, const std::vector<value_type>*);
    FactorBuffer<Factor>& connectFactorHullTRBP(const size_t&, FactorBuffer<Factor>&); 
    size_t numberOfBuffers() const;
    void propagateAll(const value_type& = 0, const bool useNormalization=false); 
    void propagateOne(const size_t&, const value_type& = 0, const bool useNormalization=false);
    void propagate(const size_t&, const value_type& = 0, const bool useNormalization=false);
    void marginal(Factor&, const bool useNormalization=true) const;
    template<class DIST> value_type distance(const size_t&) const;
    //const value_type getBound() const {return bound_;}

private:
    size_t id(const size_t&);

    Factor factor_;
    std::vector<FactorBuffer<Factor>* > outBuffer_;
    std::vector<FactorBuffer<Factor>  > inBuffer_;
    std::vector<size_t> factorIndices_;
    std::vector<size_t> variableIndices_;
    std::vector<value_type>* rho_;
    //value_type bound_;
};

/// Wrapper class for factor nodes
template<class Factor, class OP, class ACC>
class FactorHullTRBP
{
public:
    typedef Factor factor_type;
    typedef typename Factor::value_type value_type;

    FactorHullTRBP();
    size_t numberOfBuffers() const {return variableIndices_.size();}
    size_t variableIndex(size_t i) const {return variableIndices_[i];}
    void assign(const Factor*, const size_t&, std::vector<VariableHullTRBP<Factor, OP, ACC> >&, const std::vector<value_type>*);
    void propagateAll(const value_type& = 0, const bool useNormalization=true);  
    void propagateOne(const size_t&, const value_type& = 0, const bool useNormalization=true); 
    void propagate(const size_t&, const value_type& = 0, const bool useNormalization=true);
    void marginal(Factor&, const bool useNormalization=true) const;
    template<class DIST> value_type distance(const size_t&) const;
  //const value_type getBound()const {value_type b = bounds_[0]; for(size_t i=1; i<bounds_.size();++i) OP::op(bounds_[i],b); return b;}

private:
    factor_type* myFactor_;
    factor_type factor_; 
    std::vector<FactorBuffer<factor_type>* > outBuffer_;
    std::vector<FactorBuffer<factor_type> > inBuffer_;
    size_t factorIndex_;
    std::vector<value_type>* rho_;
    std::vector<size_t> variableIndices_;
  //std::vector<value_type>  bounds_;
};

template<class GM, class ACC, class DIST>
class TreeReweightedBeliefPropagation : public Inference<GM, ACC>
{
public:
    typedef GM gm_type;
    typedef ACC Accumulation;
    typedef DIST Distance;
    typedef typename gm_type::factor_type factor_type;
    typedef typename factor_type::value_type value_type;   
    typedef typename factor_type::space_type space_type;
    typedef typename space_type::state_type state_type;
    typedef typename gm_type::Operator Operator;

    struct Parameter {
        size_t maximumNumberOfSteps_;
        value_type bound_;
        value_type damping_;
        bool inferSequential_;
        std::vector<size_t> sortedNodeList_;
        std::vector<value_type> rho_; 
        bool useNormalization_;
        Parameter(const size_t& maximumNumberOfSteps = 50,
            const value_type& bound = static_cast<value_type>(0.0000001),
            const value_type& damping = static_cast<value_type>(0)
            )
            : maximumNumberOfSteps_(maximumNumberOfSteps),
            bound_(bound),
            damping_(damping),
            inferSequential_(false),
	    useNormalization_(true)
        {}
    };

    struct Message {
        size_t nodeId_;
        size_t internalMessageId_;
        Message() 
            : nodeId_(-1), internalMessageId_(-1)
        {}
        Message(const size_t& nodeId, const size_t& internalMessageId)
            : nodeId_(nodeId), internalMessageId_(internalMessageId)
        {}
    };

    // construction
    TreeReweightedBeliefPropagation(const gm_type&, const Parameter& = Parameter());

    // inference
    InferenceTermination infer(); 
    void inferAsynchronous();	
    void inferSynchronous();  
    void inferSequential();
    template<class VisitorType>
    InferenceTermination infer(VisitorType&);
    template<class VisitorType>
    void inferSynchronous(VisitorType&);
    template<class VisitorType>
    void inferAsynchronous(VisitorType&);
    template<class VisitorType>
    void inferSequential(VisitorType&);
    void propagate(const value_type& = 0); 
    InferenceTermination arg(std::vector<state_type>&, const size_t& = 1)const; 
    InferenceTermination bound(value_type&) const;

    // query
    std::string name() const;
    const gm_type& graphicalModel() const;
    InferenceTermination marginal(const size_t&, factor_type& out) const;
    InferenceTermination factorMarginal(const size_t&, factor_type& out) const;
    value_type convergenceXF() const;
    value_type convergenceFX() const;
    value_type convergence() const;
    const std::vector<value_type>& convergenceHistory() const;

private:
    const gm_type& gm_;
    Parameter parameter_;
    std::vector<FactorHullTRBP<factor_type, Operator, ACC> > factorHulls_;
    std::vector<VariableHullTRBP<factor_type, Operator, ACC> > variableHulls_;
    std::vector<value_type> convergenceHistory_;
};

// visitors

template<class TRBP>
class TreeReweightedBeliefPropagationVisitor {
public:
    typedef TRBP trbp_type;

    void operator()(const trbp_type&) const;
};

template<class TRBP>
class TreeReweightedBeliefPropagationVerboseVisitor {
public:
    typedef TRBP trbp_type;

    TreeReweightedBeliefPropagationVerboseVisitor();
    void operator()(const trbp_type&);

private:
    size_t step_;
};

// implementation of VariableHullTRBP

template<class Factor, class OP, class ACC>
VariableHullTRBP<Factor, OP, ACC>::VariableHullTRBP(){}

template<class Factor, class OP, class ACC>
inline void 
VariableHullTRBP<Factor, OP, ACC>::assign
(
    const space_type& space, 
    const size_t& variableIndex, 
    const std::vector<size_t>& factorIndices, 
    const std::vector<value_type>* rho
)
{
    rho_ = (std::vector<value_type>*)(rho);
    variableIndices_.resize(1);
    variableIndices_[0] = variableIndex;
    factorIndices_  = factorIndices;
    factor_.assign(space,variableIndices_.begin(), variableIndices_.end());
    inBuffer_.resize(factorIndices.size()); 
    outBuffer_.resize(factorIndices.size());
    // allocate input-buffer
    for(size_t j=0; j<factorIndices_.size(); ++j){
        inBuffer_[j].assign(space, variableIndices_.begin(), variableIndices_.end(), OP::template neutral<value_type>());
    } 
    //bound_ = OP::template neutral<value_type>();
}

template<class Factor, class OP, class ACC>
inline size_t 
VariableHullTRBP<Factor, OP, ACC>::id
(
    const size_t& factorIndex
)
{
    for(size_t j=0; j<numberOfBuffers(); ++j) {
        if(factorIndices_[j] == factorIndex) {
            return j;
        }
    }
    // factor index not found
    assert(NO_ARG_TEST || false); 
    return 0;
}

template<class Factor, class OP, class ACC>
inline size_t 
VariableHullTRBP<Factor, OP, ACC>::numberOfBuffers() const 
{
    return factorIndices_.size();
}

template<class Factor, class OP, class ACC>
inline FactorBuffer<Factor>& 
VariableHullTRBP<Factor, OP, ACC>::connectFactorHullTRBP
(
    const size_t& factorIndex, 
    FactorBuffer<Factor>& variableOutBuffer
)
{
    size_t j = id(factorIndex);
    outBuffer_[j] = &variableOutBuffer;
    return inBuffer_[j];
}

template<class Factor, class OP, class ACC >
void 
VariableHullTRBP<Factor, OP, ACC>::propagate
(
    const size_t& id, 
    const value_type& damping,
    const bool useNormalization
)
{
    Assert(NO_ARG_TEST || id < outBuffer_.size());

    outBuffer_[id]->toggle();
    if(factorIndices_.size() < 2) {
        // nothing to send
        return; 
    }

    // initialize neutral message
    Factor& newMessage = outBuffer_[id]->current();
    newMessage.assign(newMessage.space(), variableIndices_.begin(),
        variableIndices_.end(),OP::template neutral<value_type>());

    // operate all messages going to x, except f -> x
    for(size_t j=0; j<factorIndices_.size(); ++j) {
        if(j != id) {
            OP::wop(inBuffer_[j].current(),(*rho_)[factorIndices_[j]],newMessage);
        }
        else{
            OP::wop(inBuffer_[j].current(),(*rho_)[factorIndices_[j]]-1,newMessage);
        }
    }
 
    // damp message
    if(damping != 0) {
        Factor& oldMessage = outBuffer_[id]->old();
        OP::weightedMean(newMessage, oldMessage, 1-damping, newMessage);
    }

    /*
    value_type v;
    newMessage.template accumulate<ACC>(v);
    if(useNormalization_){
      OP::iop(v,newMessage); 
    }
    else{ 
      ;
    }
    */
    //normalize message 
    if(useNormalization)
      OP::normalize(newMessage); 
}  

template<class Factor, class OP, class ACC >
inline void 
VariableHullTRBP<Factor, OP, ACC>::propagateOne
(
    const size_t& factorIndex,
    const value_type& damping,
    const bool useNormalization
)
{
    for(size_t j=0; j<factorIndices_.size(); ++j) {
        if(factorIndices_[j] == factorIndex) {
	  propagate(j, damping, useNormalization);
        }
    }
} 

template<class Factor, class OP, class ACC>
inline void 
VariableHullTRBP<Factor, OP, ACC>::propagateAll
(
    const value_type& damping,
    const bool useNormalization
)
{
    for(size_t j=0; j<factorIndices_.size(); ++j) {
      propagate(j, damping, useNormalization);
    }
} 

template<class Factor, class OP, class ACC >
inline void 
VariableHullTRBP<Factor, OP, ACC>::marginal
(
    Factor& out,
    const bool useNormalization
) const
{
    // set out to neutral
    out.assign(factor_.space(), variableIndices_.begin(),
        variableIndices_.end(), OP::template neutral<value_type>());

    // operate all incomming messages
    for(size_t j=0; j<numberOfBuffers(); ++j) {
        OP::wop(inBuffer_[j].current(),(*rho_)[factorIndices_[j]],out); 
    } 

    // normalize output
    if(useNormalization)
      OP::normalize(out);
} 

template<class Factor,class OP, class ACC >
template<class DIST> 
inline typename Factor::value_type 
VariableHullTRBP<Factor, OP, ACC>::distance
(
    const size_t& j
) const
{
    return inBuffer_[j].template dist<DIST>();
}

// implementation of FactorHullTRBP

template<class Factor, class OP, class ACC>
FactorHullTRBP<Factor, OP, ACC>::FactorHullTRBP()
{}

template<class Factor, class OP, class ACC>
inline void 
FactorHullTRBP<Factor, OP, ACC>::assign
(
    const Factor* factor, 
    const size_t& factorIndex, 
    std::vector<VariableHullTRBP<Factor, OP, ACC> >& variableHulls,
    const std::vector<value_type>* rho
)
{
    rho_ = (std::vector<value_type>*) rho;
    myFactor_ = (Factor*)(factor);
    factor_ = *factor;
    factorIndex_ = factorIndex;

    variableIndices_.resize(factor->numberOfVariables());
    factor->variableIndices(variableIndices_.begin());

    inBuffer_.resize(factor->numberOfVariables()); 
    outBuffer_.resize(factor->numberOfVariables());

    size_t count = 0;
    for(std::vector<size_t>::iterator it = variableIndices_.begin(); it != variableIndices_.end(); ++it){
        inBuffer_[count].assign(factor->space(), it, it+1, OP::template neutral<value_type>());
        outBuffer_[count] = &(variableHulls[*it].connectFactorHullTRBP(factorIndex_, inBuffer_[count]));
        ++count;
    } 

    //bounds_.resize(factor->numberOfVariables(), OP::template neutral<value_type>());
} 

template<class Factor, class OP, class ACC >
void 
FactorHullTRBP<Factor, OP, ACC>::propagate
(
    const size_t& id, 
    const value_type& damping ,
    const bool useNormalization
)
{ 
    Assert(NO_ARG_TEST || id < outBuffer_.size());

    outBuffer_[id]->toggle();
    Factor& newMessage = outBuffer_[id]->current();
    std::vector<size_t> accVar;

    // initialize
    OP::ihop(*myFactor_,(*rho_)[factorIndex_],factor_);

    // operate
    for(size_t j=0; j<variableIndices_.size(); ++j) {
        if(j != id) {
            OP::op(inBuffer_[j].current(), factor_);
            accVar.push_back(variableIndices_[j]);
        }
    }

    // accumulation over all variables except x
    factor_.template accumulate<ACC>(accVar.begin(), accVar.end(), newMessage);

    // damp message
    if(damping != 0) {
        Factor& oldMessage = outBuffer_[id]->old();
        OP::weightedMean(newMessage, oldMessage, 1-damping, newMessage);
    }
    /*
    value_type v;
    newMessage.template accumulate<ACC>(v);
    if(useNormalization_){
      OP::iop(v,newMessage); 
    }
    else{ 
      ;
    }
    */
    // normalize message
    if(useNormalization)
      OP::normalize(newMessage);
}

template<class Factor, class OP, class ACC >
inline void 
FactorHullTRBP<Factor, OP, ACC>::propagateOne
(
    const size_t& variableIndex,
    const value_type& damping,
    const bool useNormalization
)
{
    for(size_t j=0; j<variableIndices_.size(); ++j) {
        if(variableIndices_[j] == variableIndex) {
	  propagate(j, damping,useNormalization);
        }
    }
} 

template<class Factor, class OP, class ACC >
inline void 
FactorHullTRBP<Factor, OP, ACC>::propagateAll
(
    const value_type& damping,
    const bool useNormalization
)
{
    for(size_t j=0; j<variableIndices_.size(); ++j) {
      propagate(j, damping, useNormalization);
    }
} 

template<class Factor, class OP, class ACC>
inline void 
FactorHullTRBP<Factor, OP, ACC>::marginal
(
    Factor& out,
    const bool useNormalization
) const
{
    OP::ihop(*(const_cast<Factor*>(myFactor_)),(*rho_)[factorIndex_],out);

    // include incoming messages
    for(size_t j=0; j<variableIndices_.size(); ++j){
        OP::op(inBuffer_[j].current(), out);
    }
    if(useNormalization)
      OP::normalize(out); 
} 

template<class Factor, class OP, class ACC >
template<class DIST> 
inline typename Factor::value_type 
FactorHullTRBP<Factor, OP, ACC>::distance
(
    const size_t& j
) const
{
    return  inBuffer_[j].template dist<DIST>();
}


//  Implementation of TreeReweightedBeliefPropagation  

template<class GM, class ACC, class DIST>
TreeReweightedBeliefPropagation<GM, ACC, DIST>::TreeReweightedBeliefPropagation
(
    const gm_type& gm,
    const Parameter& parameter
)
: gm_(gm), 
  parameter_(parameter)
{
    if(parameter_.sortedNodeList_.size() == 0){
        parameter_.sortedNodeList_.resize( gm.space().dimension()  );
        for(size_t i=0; i<gm.space().dimension(); ++i)
            parameter_.sortedNodeList_[i] = i;
    }
    assert(parameter_.sortedNodeList_.size() == gm.space().dimension());

    std::vector<std::vector<size_t> > variableNeighbours;

    // set rho if not set manually
    if(parameter_.rho_.size() == 0) {
        // set rho by tree decomposition
        opengm::Decomposer<gm_type> decomposer;
        std::vector<std::vector<size_t> > forestDecomposition;
        decomposer.decompose(gm_, forestDecomposition);
        parameter_.rho_.resize(gm_.numberOfFactors());
        for(size_t j=0; j<forestDecomposition.size(); ++j) {
            for(size_t k=0; k<forestDecomposition[j].size(); ++k) {
                ++parameter_.rho_[forestDecomposition[j][k]];
            }
        }
        // std::cout << "A = [";
        for(size_t j=0; j<parameter_.rho_.size(); ++j) {
            parameter_.rho_[j] /= 
                static_cast<value_type>(forestDecomposition.size());
            // std::cout << parameter_.rho_[j] << ' ';
        }
        // std::cout << "];" << std::endl;
    }
    else if(parameter_.rho_.size() != gm.numberOfFactors()) {
        throw std::runtime_error("rho not set correctly.");
    }

    // test rho
    assert( parameter_.rho_.size() == gm.numberOfFactors() );
    for(size_t i=0; i<gm.numberOfFactors(); ++i){
        if(gm_[i].space().dimension()<2){
            assert(parameter_.rho_[i] == 1); // ??? allow for numerical deviation
        }
        assert(parameter_.rho_[i]>0);
    }

    //Calculate variable neighbours
    variableNeighbours.resize(gm.space().dimension());
    for(size_t j=0; j<gm.space().dimension(); ++j) {
        variableNeighbours[j].resize(0);
    }
    for(size_t i=0; i<gm.numberOfFactors(); ++i) {
        for(size_t j=0; j<gm[i].numberOfVariables(); ++j) {
            variableNeighbours[gm[i].variableIndex(j)].push_back(i);
        }
    }

    //Set Hulls
    variableHulls_.resize(gm.space().dimension(),VariableHullTRBP<factor_type, Operator, ACC>());
    for(size_t i=0; i<gm.space().dimension(); ++i) { 
      variableHulls_[i].assign(gm.space(), i, variableNeighbours[i], &parameter_.rho_);
    }
    factorHulls_.resize(gm.numberOfFactors(),FactorHullTRBP<factor_type, Operator, ACC>());
    for(size_t i=0; i<gm.numberOfFactors(); i++){
        factorHulls_[i].assign(&(gm_[i]), i, variableHulls_, &parameter_.rho_);	
    }
}

template<class GM, class ACC, class DIST>
inline std::string 
TreeReweightedBeliefPropagation<GM, ACC, DIST>::name() const
{
    return "TRBP";
}

template<class GM, class ACC, class DIST>
inline const typename TreeReweightedBeliefPropagation<GM, ACC, DIST>::gm_type&
TreeReweightedBeliefPropagation<GM, ACC, DIST>::graphicalModel() const
{
    return gm_;
} 

template<class GM, class ACC, class DIST>
inline InferenceTermination 
TreeReweightedBeliefPropagation<GM, ACC, DIST>::infer()
{
    TreeReweightedBeliefPropagationVisitor<TreeReweightedBeliefPropagation<GM, ACC, DIST> > visitor;
    return infer(visitor);
}

template<class GM, class ACC, class DIST>
template<class VisitorType>
inline InferenceTermination 
TreeReweightedBeliefPropagation<GM, ACC, DIST>::infer
(
    VisitorType& visitor
)
{
    if(gm_.isAcyclic()) {
        inferAsynchronous(visitor);
    }
    else{
        if(parameter_.inferSequential_)
            inferSequential(visitor);
        else
            inferSynchronous(visitor);
    }
    return NORMAL;
}

/// Asynchronous message passing.
///
/// A message is sent from a variable (resp. from a factor) only if 
/// all messages to that variable (factor) have been received. 
/// Asynchronous message passing is suitable for graphical models
/// whose factor graph is acyclic.
///
template<class GM, class ACC, class DIST>
void 
TreeReweightedBeliefPropagation<GM, ACC, DIST>::inferAsynchronous()
{
    TreeReweightedBeliefPropagationVisitor<TreeReweightedBeliefPropagation<GM, ACC, DIST> > visitor;
    inferAsynchronous(visitor);
}

/// Asynchronous message passing.
///
/// A message is sent from a variable (resp. from a factor) only if 
/// all messages to that variable (factor) have been received. 
/// Asynchronous message passing is suitable for graphical models
/// whose factor graph is acyclic.
///
template<class GM, class ACC, class DIST>
template<class VisitorType>
void 
TreeReweightedBeliefPropagation<GM, ACC, DIST>::inferAsynchronous
(
    VisitorType& visitor
)
{
    assert(gm_.isAcyclic()); 
    for(size_t fac=0; fac<gm_.numberOfFactors(); ++fac) {
        assert(parameter_.rho_[fac]==1);
    }

    size_t numberOfVariables = gm_.space().dimension();
    size_t numberOfFactors = gm_.numberOfFactors();

    std::vector<std::vector<size_t> > neighboursOfVariableNodes(numberOfVariables);

    // number of messages which have not yet been recevied 
    // but are required for sending
    std::vector<std::vector<size_t> > counterVar2FacMessage(numberOfVariables);
    std::vector<std::vector<size_t> > counterFac2VarMessage(numberOfFactors);

    // list of messages which are ready to send
    std::list<Message> ready2SendVar2FacMessage;
    std::list<Message> ready2SendFac2VarMessage;

    for(size_t fac=0; fac<numberOfFactors; ++fac) {
        counterFac2VarMessage[fac].resize(gm_[fac].numberOfVariables(),gm_[fac].numberOfVariables()-1);
        for(size_t i=0; i<gm_[fac].numberOfVariables(); ++i) {
            neighboursOfVariableNodes[gm_[fac].variableIndex(i)].push_back(fac);
        }
    }
    for(size_t var=0; var<numberOfVariables; ++var) {
        counterVar2FacMessage[var].resize(neighboursOfVariableNodes[var].size());
        for(size_t i=0; i<neighboursOfVariableNodes[var].size(); ++i) {
            counterVar2FacMessage[var][i] = neighboursOfVariableNodes[var].size()-1;
        }
    }

    // find all messages which are ready for sending
    for(size_t var=0; var<numberOfVariables; ++var) {
        for(size_t i=0; i<counterVar2FacMessage[var].size(); ++i) {
            if(counterVar2FacMessage[var][i] == 0) {
                --counterVar2FacMessage[var][i];
                ready2SendVar2FacMessage.push_back(Message(var,i));
            }
        }
    }
    for(size_t fac=0; fac<numberOfFactors; ++fac) {
        for(size_t i=0; i<counterFac2VarMessage[fac].size(); ++i) {
            if(counterFac2VarMessage[fac][i] == 0) {
                --counterFac2VarMessage[fac][i];
                ready2SendFac2VarMessage.push_back(Message(fac,i));
            }
        }
    }

    // send messages
    while(ready2SendVar2FacMessage.size()>0 || ready2SendFac2VarMessage.size()>0) {
        while(ready2SendVar2FacMessage.size()>0) {
            Message m = ready2SendVar2FacMessage.front();
            size_t nodeId = m.nodeId_;
            size_t factorId = neighboursOfVariableNodes[nodeId][m.internalMessageId_];

            // send message
            variableHulls_[nodeId].propagate(m.internalMessageId_, 0, false);
            ready2SendVar2FacMessage.pop_front();

            //check if new messages can be sent
            for(size_t i=0; i<gm_[factorId].numberOfVariables(); ++i) {
                if(gm_[factorId].variableIndex(i)!=nodeId) { 
                    if(--counterFac2VarMessage[factorId][i]==0) {
                        ready2SendFac2VarMessage.push_back(Message(factorId,i));
                    }
                }
            }
        }
        while(ready2SendFac2VarMessage.size()>0) {
            Message m = ready2SendFac2VarMessage.front();
            size_t factorId = m.nodeId_;
            size_t nodeId = gm_[factorId].variableIndex(m.internalMessageId_);

            // send message
            factorHulls_[factorId].propagate(m.internalMessageId_, 0, parameter_.useNormalization_);
            ready2SendFac2VarMessage.pop_front();

            // check if new messages can be sent
            for(size_t i=0; i<neighboursOfVariableNodes[nodeId].size(); ++i){
                if(neighboursOfVariableNodes[nodeId][i]!=factorId){
                    if(--counterVar2FacMessage[nodeId][i]==0){
                        ready2SendVar2FacMessage.push_back(Message(nodeId,i));
                    }
                }
            }
        }
        visitor(*this);
    }
}

template<class GM, class ACC, class DIST>
inline void 
TreeReweightedBeliefPropagation<GM, ACC, DIST>::propagate
(
    const value_type& damping
)
{
    for(size_t i=0; i<variableHulls_.size(); ++i){
      variableHulls_[i].propagateAll(damping, false);  
    } 
    for(size_t i=0; i<factorHulls_.size(); ++i){
      factorHulls_[i].propagateAll(damping, parameter_.useNormalization_);
    } 
}

template<class GM, class ACC, class DIST>
inline void 
TreeReweightedBeliefPropagation<GM, ACC, DIST>::inferSynchronous()
{
    TreeReweightedBeliefPropagationVisitor<TreeReweightedBeliefPropagation<GM, ACC, DIST> > visitor;
    inferSynchronous(visitor);
}

template<class GM, class ACC, class DIST>
template<class VisitorType>
inline void 
TreeReweightedBeliefPropagation<GM, ACC, DIST>::inferSynchronous
(
    VisitorType& visitor
)
{
    value_type c = 0;
    value_type damping = parameter_.damping_;

    
    //Let all Factors with a order lower than 2 sending their Message
    for(size_t i=0; i<factorHulls_.size(); ++i){
        if(factorHulls_[i].numberOfBuffers()<2){
	  factorHulls_[i].propagateAll(0, parameter_.useNormalization_);
	  factorHulls_[i].propagateAll(0, parameter_.useNormalization_);//2 times to fill both buffers
        }
    }
    
    for(unsigned long n=0; n<parameter_.maximumNumberOfSteps_; ++n) {
        for(size_t i=0; i<variableHulls_.size(); ++i){
	  variableHulls_[i].propagateAll(damping, false);  
        } 
        for(size_t i=0; i<factorHulls_.size(); ++i){
	   if(factorHulls_[i].numberOfBuffers()>=2)//Messages from factors of the order lower than 2 does not change!
                factorHulls_[i].propagateAll(damping, parameter_.useNormalization_);
        } 
        visitor(*this);
        c = convergence();
        convergenceHistory_.push_back(c);
        if(c < parameter_.bound_) {
            break;
        }
    }
}

// Sequential Messages according to Vladimir Kolmogorov (TRW-S) and Tappen (BP-S)
// NOTE: This algorithms are designed for models of order 2, we are not sure if our implementation
//       hold the convergence properties for higher order graphs
template<class GM, class ACC, class DIST>
inline void 
TreeReweightedBeliefPropagation<GM, ACC, DIST>::inferSequential()
{
    TreeReweightedBeliefPropagationVisitor<TreeReweightedBeliefPropagation<GM, ACC, DIST> > visitor;
    inferSequential(visitor);
}

// Sequential Messages according to Vladimir Kolmogorov (TRW-S) and Tappen (BP-S)
// NOTE: This algorithms are designed for models of order 2, we are not sure if our implementation
//       hold the convergence properties for higher order graphs
template<class GM, class ACC, class DIST>
template<class VisitorType>
inline void 
TreeReweightedBeliefPropagation<GM, ACC, DIST>::inferSequential
(
    VisitorType& visitor
)
{
    assert(parameter_.sortedNodeList_.size() == gm_.space().dimension());
    value_type damping = parameter_.damping_;

    //set nodeOrder
    std::vector<size_t> nodeOrder( gm_.space().dimension() );
    for(size_t o=0; o<gm_.space().dimension(); ++o){ 
        nodeOrder[parameter_.sortedNodeList_[o]]=o;
    }
    

    //Let all Factors with a order lower than 2 sending their Message
    for(size_t i=0; i<factorHulls_.size(); ++i){
        if(factorHulls_[i].numberOfBuffers()<2){
            factorHulls_[i].propagateAll(0, parameter_.useNormalization_); 
            factorHulls_[i].propagateAll(0, parameter_.useNormalization_);//2 times to fill both buffers
        }
    }
    

    //The following Code is not optimized and maybe too slow for small factors 
    for(unsigned long n=0; n<parameter_.maximumNumberOfSteps_; ++n){
        //in increasing ordering 
        for(size_t o=0; o<gm_.space().dimension(); ++o){
            size_t node = parameter_.sortedNodeList_[o]; 

            //Update messages from Variable
            variableHulls_[node].propagateAll(damping,false);

            //Send all messages form factor to variable with a nodeorder > o
            //which require only messages from nodes with a nodeorder <= o 
            //and are not passed in this run (one has nodeorder = o
            for(size_t i=0; i<factorHulls_.size(); ++i){
	      if(factorHulls_[i].numberOfBuffers()>=2){
                    bool sendMessage = true;
                    bool setInternalId = false;
                    // size_t internalId;
                    for(size_t j=0; j<factorHulls_[i].numberOfBuffers(); ++j){
                        size_t factorNode = factorHulls_[i].variableIndex(j);
                        if(factorNode==node){
                            factorHulls_[i].propagateAll(damping, parameter_.useNormalization_);
                        }
                    }
	      }
            }

        }
        //in decreasing ordering 
        for(size_t oo=0; oo<gm_.space().dimension(); ++oo){
            size_t o    = gm_.space().dimension()-1-oo;
            size_t node = parameter_.sortedNodeList_[o];

            //Update messages from Variable 
            variableHulls_[node].propagateAll(damping,false);

            //Send all messages form factor to variable U with a nodeorder[X] < o
            //which require one messages from node V with a nodeorder[V] = o
            //and no message from a node W with  nodeorder[X]<nodeorder[W]<o
            for(size_t i=0; i<factorHulls_.size(); ++i){
	      if(factorHulls_[i].numberOfBuffers()>=2){
                    bool sendMessage = false;
                    bool setInternalId = true;
                    // size_t internalId;
                    for(size_t j=0; j<factorHulls_[i].numberOfBuffers(); ++j){
                        size_t factorNode = factorHulls_[i].variableIndex(j);
                        if(factorNode==node)
                            factorHulls_[i].propagateAll(damping, parameter_.useNormalization_);
                    }
	      }
            }
        }
        visitor(*this);
    }  
}

template<class GM, class ACC, class DIST>
inline InferenceTermination 
TreeReweightedBeliefPropagation<GM, ACC, DIST>::marginal
(
    const size_t& variableIndex,
    factor_type& out
) const
{
    Assert(NO_ARG_TEST || variableIndex < variableHulls_.size());
    variableHulls_[variableIndex].marginal(out, parameter_.useNormalization_);
    return NORMAL;
}

template<class GM, class ACC, class DIST>
inline InferenceTermination 
TreeReweightedBeliefPropagation<GM, ACC, DIST>::factorMarginal
(
    const size_t& factorIndex,
    factor_type &out
) const
{
    Assert(NO_ARG_TEST || factorIndex < factorHulls_.size());
    factorHulls_[factorIndex].marginal(out, parameter_.useNormalization_);
    return NORMAL;
}

// maximum distance among all pairs of variable to factor messages
// between the previous and the current belief propagation step
template<class GM, class ACC, class DIST>
inline typename TreeReweightedBeliefPropagation<GM, ACC, DIST>::value_type
TreeReweightedBeliefPropagation<GM, ACC, DIST>::convergenceXF() const
{ 
    value_type result = 0;
    for(size_t j=0; j<factorHulls_.size(); ++j) {
        for(size_t i=0; i<factorHulls_[j].numberOfBuffers(); ++i) {
            value_type d = factorHulls_[j].template distance<DIST>(i);
            if(d > result) {
                result = d;
            }
        }
    }
    return result;
}

// maximum distance among all pairs of factor to variable messages
// between the previous and the current belief propagation step
template<class GM, class ACC, class DIST>
inline typename TreeReweightedBeliefPropagation<GM, ACC, DIST>::value_type
TreeReweightedBeliefPropagation<GM, ACC, DIST>::convergenceFX() const
{ 
    value_type result = 0;
    for(size_t j=0; j<variableHulls_.size(); ++j) {
        for(size_t i=0; i<variableHulls_[j].numberOfBuffers(); ++i) {
            value_type d = variableHulls_[j].template distance<DIST>(i);
            if(d > result) {
                result = d;
            }
        }
    }
    return result;
}

template<class GM, class ACC, class DIST>
inline typename TreeReweightedBeliefPropagation<GM, ACC, DIST>::value_type
TreeReweightedBeliefPropagation<GM, ACC, DIST>::convergence() const
{
    return convergenceXF();
}

template<class GM,class ACC,class DIST >
inline InferenceTermination 
TreeReweightedBeliefPropagation<GM, ACC, DIST>::arg
(
    std::vector<state_type>& conf, 
    const size_t& N
)const
{
    if(N != 1) {
        // belief propagation cannot return k-th optimal configuration
        Assert(0 == 1);
        return UNKNOWN;
    }
    else { 
      return modeFromMarginal(conf);
    }
}


template<class GM,class ACC,class DIST >
inline const std::vector<typename TreeReweightedBeliefPropagation<GM, ACC, DIST>::value_type>&
TreeReweightedBeliefPropagation<GM, ACC, DIST>::convergenceHistory() const
{
    return convergenceHistory_;
}

template<class GM,class ACC,class DIST > 
InferenceTermination
TreeReweightedBeliefPropagation<GM, ACC, DIST>::bound(value_type& bound) const
{
  bound = ACC::template ineutral<value_type>();
  return UNKNOWN;  
}

// implementation of visitors

template<class TRBP>
inline void
TreeReweightedBeliefPropagationVisitor<TRBP>::operator()
(
    const typename TreeReweightedBeliefPropagationVisitor<TRBP>::trbp_type& trbp
) const
{}

template<class TRBP>
inline 
TreeReweightedBeliefPropagationVerboseVisitor<TRBP>::TreeReweightedBeliefPropagationVerboseVisitor()
: step_(0)
{}

template<class TRBP>
inline void
TreeReweightedBeliefPropagationVerboseVisitor<TRBP>::operator()
(
    const typename TreeReweightedBeliefPropagationVerboseVisitor<TRBP>::trbp_type& trbp
)
{
    ++step_;
    std::vector<size_t> state;
    std::vector<size_t> knownVariables; 
    std::vector<size_t> knownStates; 
    trbp.constrainedOptimum(knownVariables, knownStates, state);
    std::cout << "step " << step_ 
        << ": E=" << trbp.graphicalModel().evaluate(state)
        << ", c=" << trbp.convergence()
	<< std::endl;
}

} // namespace opengm

#endif // #ifndef OPENGM_BELIEFPROPAGATION_HXX
