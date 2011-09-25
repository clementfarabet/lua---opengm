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
#ifndef OPENGM_GRAPHICALMODEL_HXX
#define OPENGM_GRAPHICALMODEL_HXX

#include <vector>
#include <set>
#include <queue>

#include "opengm/opengm.hxx"

namespace opengm {

template<class Factor, class OP>
class GraphicalModel {
public:
    typedef Factor factor_type;   
    typedef typename factor_type::value_type value_type;
    typedef typename factor_type::space_type space_type;
    typedef typename space_type::state_type state_type;
    typedef OP Operator;

    // construction
    GraphicalModel(const size_t& = 0);
    template<class FactorIterator>
        GraphicalModel(FactorIterator, FactorIterator); // copies factors 

    // query
    const space_type& space() const;
    size_t numberOfFactors() const;
    size_t factorOrder() const;
    bool isAcyclic() const;

    // access to factors
    Factor& operator[](const size_t&);
    const Factor& operator[](const size_t&) const;

    // manipulation
    void clear();
    size_t addFactor(const Factor&);
    template<class Iterator>
        value_type evaluate(Iterator) const;
    template<class IndexIterator, class ValueIterator> 
        void introduceEvidence(IndexIterator, IndexIterator, ValueIterator);

private:
    std::vector<Factor> factors_;
};

template<class Factor, class OP>
GraphicalModel<Factor, OP>::GraphicalModel
(
    const size_t& numberOfFactors
)
:   factors_(std::vector<Factor>())
{
    if(numberOfFactors != 0) {
        factors_.reserve(numberOfFactors);
    }
}

template<class Factor, class OP>
template<class FactorIterator>
GraphicalModel<Factor, OP>::GraphicalModel
(
    FactorIterator begin,
    FactorIterator end
) : factors_(std::vector<Factor>())
{
    ptrdiff_t s = std::distance(begin, end);
    if(!NO_ARG_TEST) {
        Assert(s >= 0);
        FactorIterator tmp = begin;
        while(tmp != end) {
            Assert(NO_ARG_TEST || &(tmp->space()) == &(begin->space()) );
            ++tmp;
        }
    }
    factors_.resize(s);
    std::copy(begin, end, factors_.begin());
}

template<class Factor, class OP>
inline const typename GraphicalModel<Factor, OP>::space_type& 
GraphicalModel<Factor, OP>::space() const
{
    Assert(NO_DEBUG || numberOfFactors() != 0);
    return factors_[0].space();
}

template<class Factor, class OP>
inline size_t 
GraphicalModel<Factor, OP>::numberOfFactors() const
{
    return factors_.size();
}

template<class Factor, class OP>
inline const Factor& 
GraphicalModel<Factor, OP>::operator[]
(
    const size_t& index
    ) const
{
    Assert(NO_ARG_TEST || index < factors_.size());
    return factors_[index];
}

template<class Factor, class OP>
inline Factor& 
GraphicalModel<Factor, OP>::operator[]
(
    const size_t& index
    )
{
    Assert(NO_ARG_TEST || index < factors_.size());
    return factors_[index];
}

template<class Factor, class OP>
inline size_t 
GraphicalModel<Factor, OP>::addFactor
(
    const Factor& factor
)
{
    if(factors_.size() == 0) {
        factors_.push_back(factor);
        return 0;
    }
    else {
        Assert(NO_ARG_TEST || &space() == &factor.space());
        factors_.push_back(factor); // copy
        return numberOfFactors()-1; // index of added factor
    }
}

template<class Factor, class OP>
inline void 
GraphicalModel<Factor, OP>::clear()
{
    factors_.clear();
}

template<class Factor, class OP>
template<class IndexIterator, class ValueIterator>
inline void 
GraphicalModel<Factor, OP>::introduceEvidence
(
    IndexIterator stateBegin,
    IndexIterator stateEnd,
    ValueIterator value
)
{
    for(size_t j=0; j<numberOfFactors(); ++j) {
        factors_[j].fixVariables(stateBegin, stateEnd, value);
    }
}

template<class Factor, class OP>
bool GraphicalModel<Factor, OP>::isAcyclic() const
{
    if(factors_.size() == 0) {
        return true;
    }

    // construct variable adjacency graph
    std::vector<std::set<size_t> > adjacency(space().dimension() + numberOfFactors());
    for(size_t f=0; f<numberOfFactors(); ++f) {
        for(size_t v=0; v<factors_[f].numberOfVariables(); ++v) {
            size_t p = space().dimension() + f; // factor
            size_t q = factors_[f].variableIndex(v);
            adjacency[p].insert(q);
            adjacency[q].insert(p);
        }
    }

    // search for loops
    std::vector<size_t> parent(space().dimension() + numberOfFactors(), -1);
    std::queue<size_t> queue;
    for(size_t v=0; v<space().dimension(); ++v) {
        if(parent[v] != -1) {
            continue;
        }
        else {
            parent[v] = v;
            queue.push(v);
            // std::cout << v << " -> " << v << std::endl;
            while(!queue.empty()) {
                size_t w = queue.front();
                queue.pop();
                for(std::set<size_t>::const_iterator it = adjacency[w].begin();
                it != adjacency[w].end(); ++it) {
                    if(parent[*it] == -1) {
                        parent[*it] = w;
                        queue.push(*it);
                        // std::cout << w << " -> " << *it << std::endl;
                    }
                    else {
                        if(parent[w] != *it) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
}

template<class Factor, class OP>
inline size_t 
GraphicalModel<Factor, OP>::factorOrder() const
{
    size_t factorOrder = 0;
    for(size_t i=0; i<numberOfFactors(); i++) {
        if(factors_[i].numberOfVariables() > factorOrder)
            factorOrder = factors_[i].numberOfVariables();
    } 
    return factorOrder;
}

template<class Factor, class OP>
template<class Iterator>
typename GraphicalModel<Factor, OP>::value_type
GraphicalModel<Factor, OP>::evaluate
(
    Iterator it
) const
{
    value_type v;
    Operator::neutral(v);
    for(size_t j=0; j<numberOfFactors(); ++j) {
        std::vector<size_t> factor_state(factors_[j].numberOfVariables(), static_cast<size_t>(0));
        for(size_t i=0; i<factors_[j].numberOfVariables(); ++i) {
            factor_state[i] = it[factors_[j].variableIndex(i)];
        }
        Operator::op(factors_[j](factor_state.begin()), v);
    }  
    return v;
}

} // namespace opengm

#endif // #ifndef OPENGM_GRAPHICALMODEL_HXX
