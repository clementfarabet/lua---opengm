/// OpenGM. Copyright (c) 2010 by Bjoern Andres and Joerg Hendrik Kappes.
///
/// This software was developed by Bjoern Andres and Joerg Hendrik Kappes.
/// Enquiries shall be directed to:
/// bjoern.andres@iwr.uni-heidelberg.de, kappes@math.uni-heidelberg.de
///
/// Author(s) of this file: Bjoern Andres, Joerg Hendrik Kappes
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
#ifndef OPENGM_INFERENCE_HXX
#define OPENGM_INFERENCE_HXX

#include <vector>
#include <string>
#include <list>

#include "opengm/opengm.hxx"

namespace opengm {  

enum InferenceTermination { 
    UNKNOWN=0, NORMAL=1, TIMEOUT=2, CONVERGENCE=3 
};

template <class GM, class ACCUMULATOR>
class Inference
{
public:
    typedef GM gm_type;
    typedef ACCUMULATOR Accumulator;
    typedef typename gm_type::factor_type factor_type;
    typedef typename factor_type::space_type space_type; 
    typedef typename factor_type::value_type value_type;
    typedef typename space_type::state_type state_type;

    // abstract member functions
    virtual std::string name() const = 0;
    virtual const GM& graphicalModel() const = 0;
    virtual InferenceTermination infer() = 0;
    virtual InferenceTermination arg(std::vector<state_type>&, const size_t& = 1) const = 0;

    // member functions with default definition
    virtual InferenceTermination args(std::vector<std::vector<state_type> >&) const;
    virtual InferenceTermination marginal(const size_t&, factor_type&) const;
    virtual InferenceTermination factorMarginal(const size_t&, factor_type&) const;
    virtual InferenceTermination bound(value_type&) const;

    InferenceTermination constrainedOptimum(std::vector<size_t>&, 
        std::vector<state_type>&, std::vector<state_type>&) const;
    InferenceTermination modeFromMarginal(std::vector<state_type>&) const;
    InferenceTermination modeFromFactorMarginal(std::vector<state_type>&) const;
};

template<class GM, class ACCUMULATOR>
inline InferenceTermination 
Inference<GM, ACCUMULATOR>::args
(
    std::vector<std::vector<state_type> >& out
)const
{
    return UNKNOWN;
}

template<class GM, class ACCUMULATOR>
inline InferenceTermination 
Inference<GM, ACCUMULATOR>::marginal
(
    const size_t& variableIndex, 
    factor_type& out
) const
{
    return UNKNOWN;
}

template<class GM, class ACCUMULATOR>
inline InferenceTermination 
    Inference<GM, ACCUMULATOR>::factorMarginal
    (
    const size_t& factorIndex, 
    factor_type& out
    ) const
{
    return UNKNOWN;
}

template<class GM, class ACCUMULATOR>
inline InferenceTermination 
Inference<GM, ACCUMULATOR>::bound
(
    value_type& out
) const
{
    return UNKNOWN;
}

template<class GM, class ACCUMULATOR>
InferenceTermination 
Inference<GM, ACCUMULATOR>::constrainedOptimum
(
    std::vector<size_t>& knownVariables, 
    std::vector<state_type>& knownStates, 
    std::vector<state_type>& conf
) const
{ 
    const GM& gm = graphicalModel();
    const space_type& space = gm.space();
    size_t numberOfNodes = space.dimension(); 
    size_t numberOfUnknowns = numberOfNodes;
    std::vector<size_t> unknown(numberOfNodes, true);   
    conf.resize(numberOfNodes);

    // set known variables
    for(size_t i=0; i<knownStates.size(); ++i){
        conf[knownVariables[i]] = knownStates[i];
        unknown[knownVariables[i]] = false;
        --numberOfUnknowns;
    }

    // make list of Factors
    std::list<size_t> factorList;
    for(size_t i=0; i<gm.numberOfFactors(); ++i) {
        const factor_type& f = gm[i];
        size_t nvar = f.numberOfVariables();
        if(nvar >= 2) {
            size_t numUnknowns = 0;
            size_t numKnowns  = 0;
            for(size_t c=0; c<nvar; ++c) {
                if(unknown[f.variableIndex(c)]) {
                    ++numUnknowns;
                }
                else {
                    ++numKnowns;
                }
            }
            if(numUnknowns > 0 && numKnowns+numUnknowns >= 2) {
                factorList.push_back(i);
            }
        }
    }

    // find a factor which connects unknown and known variables
    std::list<size_t>::iterator it = factorList.begin();
    bool foundFactor;
    size_t foundFactorId;
    while(numberOfUnknowns > 0) {
        foundFactor = false;
        if(factorList.size() > 0) {
            if(it == factorList.end()) {
                it = factorList.begin();
            }
            assert(it != factorList.end());
            size_t end = (*it);
            do {
                size_t numUnknowns = 0;
                size_t numKnowns  = 0;
                const factor_type& f = gm[(*it)];
                for(size_t c=0; c<f.numberOfVariables(); ++c) {
                    if(unknown[f.variableIndex(c)]) {
                        ++numUnknowns;
                    }
                    else {
                        ++numKnowns;
                    }
                }
                if(numUnknowns > 0 && numKnowns > 0) {
                    foundFactor = true;
                    foundFactorId = (*it);
                    it = factorList.erase(it);
                    break;
                }
                else if(numUnknowns == 0) {
                    if(*it==end){
                        it = factorList.erase(it);  
                        break;
                    }
                    else{
                        it = factorList.erase(it);
                    }
                    if(factorList.size() == 0) {
                        break;
                    }
                }
                else {
                    ++it;
                }
                if(it == factorList.end()) {
                    it=factorList.begin();
                }
            } while((*it) != end);
        }

        if(foundFactor) {
            std::vector<size_t> accVariables;
            std::vector<size_t> accStates;
            std::vector<size_t> unknownVariables;	
            const factor_type& f = gm[foundFactorId];
            for(size_t c=0; c<f.numberOfVariables(); ++c) {
                if(!unknown[f.variableIndex(c)]) {
                    accVariables.push_back(f.variableIndex(c));
                    accStates.push_back(conf[f.variableIndex(c)]);
                }
                else {
                    unknownVariables.push_back(f.variableIndex(c));
                }
            } 
            value_type value;
            std::vector<size_t> state(unknownVariables.size()); 
            factor_type beliefX;
            if(NORMAL != factorMarginal(foundFactorId, beliefX))
                return UNKNOWN;
            beliefX.fixVariables(accVariables.begin(), accVariables.end(), accStates.begin());
            beliefX.template accumulate<Accumulator>(value,state);
            for(size_t i=0; i<state.size(); ++i) {
                conf[unknownVariables[i]] = state[i];
                unknown[unknownVariables[i]] = false;
                --numberOfUnknowns;
            }
        }
        else { 
            for(foundFactorId=0; foundFactorId<numberOfNodes; ++foundFactorId) {
                if(unknown[foundFactorId]) break;
            }
            value_type value;
            std::vector<size_t> state(1); 
            factor_type beliefX; 
            if(NORMAL != marginal(foundFactorId, beliefX))
                return UNKNOWN;
            beliefX.template accumulate<Accumulator>(value,state);
            conf[foundFactorId] = state[0];
            unknown[foundFactorId] = false;
            --numberOfUnknowns;
        }
    } 
    return NORMAL;
}

template<class GM, class ACCUMULATOR>
InferenceTermination 
Inference<GM, ACCUMULATOR>::modeFromMarginal
( 
    std::vector<state_type>& conf
) const
{ 
    const GM&         gm = graphicalModel();
    const space_type& space = gm.space();
    size_t            numberOfNodes = space.dimension(); 

    conf.resize(gm.space().dimension());
    factor_type out;

    for(size_t node=0; node<numberOfNodes; ++node) {
        InferenceTermination term = marginal(node, out);
        if(NORMAL != term) {
            return term; 
        }
        value_type value = out(0);
        size_t state = 0;
        for(size_t i=1; i<gm.space().numberOfStates(node); ++i) {
            if(ACCUMULATOR::bop(out(i), value)) {
                value = out(i);
                state = i;
            }
        }
        conf[node] = state;
    }
    return NORMAL;  
} 

template<class GM, class ACCUMULATOR>
InferenceTermination 
Inference<GM, ACCUMULATOR>::modeFromFactorMarginal
( 
    std::vector<state_type>& conf
) const
{ 
    const GM&         gm = graphicalModel();
    const space_type& space = gm.space();
    size_t            numberOfNodes = space.dimension(); 

    conf.resize(gm.space().dimension());
    std::vector<size_t> none;
    none.resize(0);
    constrainedOptimum(none, none, conf);
    return NORMAL;  
}

} // namespace opengm

#endif // #ifndef OPENGM_INFERENCE_HXX
