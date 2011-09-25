/// OpenGM. Copyright (c) 2010 by Bjoern Andres and Joerg Hendrik Kappes.
///
/// This software was developed by Bjoern Andres and Joerg Hendrik Kappes.
/// Enquiries shall be directed to:
/// bjoern.andres@iwr.uni-heidelberg.de, kappes@math.uni-heidelberg.de
///
/// Author(s) of this file: Bjoern Andres
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
#ifndef OPENGM_MOVEMAKER_HXX
#define OPENGM_MOVEMAKER_HXX

#include <algorithm>

namespace opengm {

// move maker inference algorithm
template<class GM>
class Movemaker {
public:
    typedef GM gm_type;
    typedef typename gm_type::factor_type factor_type;
    typedef typename factor_type::value_type value_type;   
    typedef typename factor_type::space_type space_type;
    typedef typename space_type::state_type state_type;
    typedef typename gm_type::Operator Operator;
    typedef typename std::vector<state_type>::const_iterator state_iterator;

    Movemaker(const gm_type&); // initial state: all zeros
    template<class StateIterator>
        Movemaker(const gm_type&, StateIterator); // sets initial state

    // query
    value_type energy() const;
    template<class IndexIterator, class StateIterator>
        value_type energyAfterMove(IndexIterator, IndexIterator, StateIterator);
    const state_type& state(const size_t&) const;
    state_iterator stateBegin() const;
    state_iterator stateEnd() const;

    // manipulation
    template<class IndexIterator, class StateIterator>
        value_type move(IndexIterator, IndexIterator, StateIterator);

private:
    const gm_type& gm_;
    std::vector<std::set<size_t> > factorsOfVariable_;
    std::vector<state_type> state_;
    std::vector<state_type> stateBuffer_; // always equal to state_ (invariant)
    value_type energy_; // energy of state state_ (invariant)
};

// implementation

template<class GM>
Movemaker<GM>::Movemaker
(
    const typename Movemaker<GM>::gm_type& gm
) 
:   gm_(gm), 
    factorsOfVariable_(gm.space().dimension()),
    state_(gm.space().dimension()),
    stateBuffer_(gm.space().dimension()),
    energy_(gm.evaluate(state_.begin()))
{
    for(size_t f=0; f<gm.numberOfFactors(); ++f) {
        for(size_t v=0; v<gm[f].numberOfVariables(); ++v) {
            factorsOfVariable_[gm[f].variableIndex(v)].insert(f);
        }
    }
}

// sets initial state
template<class GM>
template<class StateIterator>
Movemaker<GM>::Movemaker
(
    const typename Movemaker<GM>::gm_type& gm,
    StateIterator it
)
:   gm_(gm), 
    factorsOfVariable_(gm.space().dimension()),
    state_(gm.space().dimension()),
    stateBuffer_(gm.space().dimension()),
    energy_(gm.evaluate(it)) // fails if *it is out of bounds
{
    for(size_t j=0; j<gm.space().dimension(); ++j, ++it) {
        state_[j] = *it;
        stateBuffer_[j] = *it;
    }
    for(size_t f=0; f<gm.numberOfFactors(); ++f) {
        for(size_t v=0; v<gm[f].numberOfVariables(); ++v) {
            factorsOfVariable_[gm[f].variableIndex(v)].insert(f);
        }
    }
}

template<class GM>
inline typename Movemaker<GM>::value_type
Movemaker<GM>::energy() const
{
    return energy_;
}

template<class GM>
template<class IndexIterator, class StateIterator>
typename Movemaker<GM>::value_type 
Movemaker<GM>::energyAfterMove
(
    IndexIterator begin, 
    IndexIterator end, 
    StateIterator destinationState
)
{
    if(!NO_ARG_TEST) {
        gm_.space().isValidIndexSequence(begin, end);
    }

    // set stateBuffer_ to destinationState, and determine factors to recompute
    std::set<size_t> factorsToRecompute;
    for(IndexIterator it = begin; it != end; ++it, ++destinationState) {
        if(state_[*it] != *destinationState) {
            stateBuffer_[*it] = *destinationState;
            std::set<size_t> tmpSet;
            std::set_union(factorsToRecompute.begin(), factorsToRecompute.end(),
                factorsOfVariable_[*it].begin(), factorsOfVariable_[*it].end(),
                std::inserter(tmpSet, tmpSet.begin()) );
            factorsToRecompute.swap(tmpSet);
        }
    }

    // ??? consider buffering the values of ALL factors at the current state!
    value_type destinationValue = energy_;
    for(std::set<size_t>::const_iterator it = factorsToRecompute.begin(); it != factorsToRecompute.end(); ++it) {
        // std::cout << "updating factor " << *it << std::endl;
        // determine current and destination state of the current factor
        std::vector<size_t> currentFactorState(gm_[*it].numberOfVariables());
        std::vector<size_t> destinationFactorState(gm_[*it].numberOfVariables());
        for(size_t j=0; j<gm_[*it].numberOfVariables(); ++j) {
            currentFactorState[j] = state_[gm_[*it].variableIndex(j)];
            destinationFactorState[j] = stateBuffer_[gm_[*it].variableIndex(j)];
        }
        // update destinationValue
        // ??? replace by operation
        destinationValue += gm_[*it](destinationFactorState.begin());
        destinationValue -= gm_[*it](currentFactorState.begin());
    }

    // restore stateBuffer_
    for(IndexIterator it = begin; it != end; ++it) {
        stateBuffer_[*it] = state_[*it];
    }

    return destinationValue;
}

template<class GM>
template<class IndexIterator, class StateIterator>
inline typename Movemaker<GM>::value_type 
Movemaker<GM>::move
(
    IndexIterator begin, 
    IndexIterator end, 
    StateIterator sit
)
{
    energy_ = energyAfterMove(begin, end, sit); // tests for assertions
    while(begin != end) {
        state_[*begin] = *sit;
        stateBuffer_[*begin] = *sit;
        ++begin;
        ++sit;
    }
    return energy_;
}

template<class GM>
inline const typename Movemaker<GM>::state_type& 
Movemaker<GM>::state
(
    const size_t& variableIndex
) const
{
    Assert(NO_ARG_TEST || variableIndex < state_.size());
    return state_[variableIndex];
}

template<class GM>
inline typename Movemaker<GM>::state_iterator 
Movemaker<GM>::stateBegin() const
{
    return state_.begin();
}

template<class GM>
inline typename Movemaker<GM>::state_iterator 
Movemaker<GM>::stateEnd() const
{
    return state_.end();
}

} // namespace opengm

#endif // #ifndef OPENGM_MOVEMAKER_HXX
