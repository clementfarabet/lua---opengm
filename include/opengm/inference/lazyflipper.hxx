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
#ifndef OPENGM_LAZYFLIPPER_HXX
#define OPENGM_LAZYFLIPPER_HXX

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <stdexcept>

#include "opengm/minimizer.hxx"
#include "opengm/inference/inference.hxx"
#include "opengm/inference/movemaker.hxx"
#include "opengm/inference/lazyflipper_tagging.hxx"
#include "opengm/inference/lazyflipper_adjacency.hxx"
#include "opengm/inference/lazyflipper_forest.hxx"

namespace opengm {

template<class GM, class ACCUMULATOR = opengm::Minimizer>
class LazyFlipper : Inference<GM, ACCUMULATOR> {
public:
    typedef GM gm_type;
    typedef ACCUMULATOR Accumulator;
    typedef typename gm_type::factor_type factor_type;
    typedef typename factor_type::value_type value_type;   
    typedef typename factor_type::space_type space_type;
    typedef typename space_type::state_type state_type;
    typedef typename gm_type::Operator Operator;
    typedef size_t VariableIndex;
    typedef Forest<VariableIndex> SubgraphForest;
    typedef size_t SubgraphForestNode;
    static const SubgraphForestNode NONODE = SubgraphForest::NONODE;

    // construction
    LazyFlipper(const gm_type&, const size_t& = 2);
    template<class StateIterator>
        LazyFlipper(const gm_type&, const size_t&, StateIterator);

    // query
    std::string name() const;
    const gm_type& graphicalModel() const;
    const size_t& maxSubgraphSize() const;
    value_type optimum() const;
    
    // manipulation
    void setMaxSubgraphSize(const size_t&);

    // inference
    InferenceTermination infer();
    template<class VisitorType>
        InferenceTermination infer(VisitorType&);
    InferenceTermination arg(std::vector<state_type>&, const size_t& = 1)const;

private:
    SubgraphForestNode appendVariableToPath(SubgraphForestNode);
    SubgraphForestNode generateFirstPathOfLength(const size_t&);
    SubgraphForestNode generateNextPathOfSameLength(SubgraphForestNode);
    void activateInfluencedVariables(SubgraphForestNode, const size_t&);
    void deactivateAllVariables(const size_t&);
    SubgraphForestNode firstActivePath(const size_t&);
    SubgraphForestNode nextActivePath(SubgraphForestNode, const size_t&);
    value_type energyAfterFlip(SubgraphForestNode);
    void flip(SubgraphForestNode);

    const gm_type& gm_;
    Adjacency variableAdjacency_;
    Movemaker<gm_type> movemaker_;
    Tagging<bool> activation_[2];
    SubgraphForest subgraphForest_;
    size_t maxSubgraphSize_;
};

// visitors

template<class LAZYFLIPPER>
class LazyFlipperVisitor {
public:
    typedef LAZYFLIPPER lazy_flipper_type;
    typedef typename lazy_flipper_type::value_type value_type;

    template<class StateIterator>
        void operator()(const lazy_flipper_type&, 
            StateIterator, StateIterator, 
            const value_type&, const size_t&, const size_t&) const;
};

template<class LAZYFLIPPER>
class LazyFlipperVerboseVisitor {
public:
    typedef LAZYFLIPPER lazy_flipper_type;
    typedef typename lazy_flipper_type::value_type value_type;

    LazyFlipperVerboseVisitor();
    template<class StateIterator>
        void operator()(const lazy_flipper_type&, 
            StateIterator, StateIterator, 
            const value_type&, const size_t&, const size_t&);

private:
    size_t step_;
};

// implementation of LazyFlipper

template<class GM, class ACCUMULATOR>
inline 
LazyFlipper<GM, ACCUMULATOR>::LazyFlipper
(
    const gm_type& gm,
    const size_t& maxSubgraphSize
)
:   gm_(gm), 
    variableAdjacency_(Adjacency(gm.space().dimension())),
    movemaker_(Movemaker<GM>(gm)),
    subgraphForest_(SubgraphForest()),
    maxSubgraphSize_(2)
{
    if(gm_.space().dimension() == 0) {
        throw std::runtime_error("The graphical model has no variables.");
    }
    for(size_t j=0; j<gm_.space().dimension(); ++j) {
        if(gm_.space().numberOfStates(j) != 2) {
            throw std::runtime_error("All variables must be binary.");
        }
    }
    setMaxSubgraphSize(maxSubgraphSize);

    // initialize activation_
    activation_[0].append(gm.space().dimension());
    activation_[1].append(gm.space().dimension());

    // initialize variableAdjacency_
    for(size_t j=0; j<gm_.numberOfFactors(); ++j) {
        const factor_type& factor = gm_[j];
        for(size_t m=0; m<factor.numberOfVariables(); ++m) {
            for(size_t n=m+1; n<factor.numberOfVariables(); ++n) {
                variableAdjacency_.connect(factor.variableIndex(m), factor.variableIndex(n));
            }
        }
    }
}

// ??? get rid of redundancy with other constructor
template<class GM, class ACCUMULATOR>
template<class StateIterator>
inline 
LazyFlipper<GM, ACCUMULATOR>::LazyFlipper
(
    const gm_type& gm,
    const size_t& maxSubgraphSize,
    StateIterator it
) 
:   gm_(gm), 
    variableAdjacency_(Adjacency(gm.space().dimension())),
    movemaker_(Movemaker<GM>(gm, it)),
    subgraphForest_(SubgraphForest()),
    maxSubgraphSize_(2)
{
    if(gm_.space().dimension() == 0) {
        throw std::runtime_error("The graphical model has no variables.");
    }
    for(size_t j=0; j<gm_.space().dimension(); ++j) {
        if(gm_.space().numberOfStates(j) != 2) {
            throw std::runtime_error("All variables must be binary.");
        }
    }
    setMaxSubgraphSize(maxSubgraphSize);

    // initialize activation_
    activation_[0].append(gm.space().dimension());
    activation_[1].append(gm.space().dimension());

    // initialize variableAdjacency_
    for(size_t j=0; j<gm_.numberOfFactors(); ++j) {
        const factor_type& factor = gm_[j];
        for(size_t m=0; m<factor.numberOfVariables(); ++m) {
            for(size_t n=m+1; n<factor.numberOfVariables(); ++n) {
                variableAdjacency_.connect(factor.variableIndex(m), factor.variableIndex(n));
            }
        }
    }
}

template<class GM, class ACCUMULATOR>
inline std::string 
LazyFlipper<GM, ACCUMULATOR>::name() const
{
    return "LazyFlipper";
}

template<class GM, class ACCUMULATOR>
inline const typename LazyFlipper<GM, ACCUMULATOR>::gm_type& 
LazyFlipper<GM, ACCUMULATOR>::graphicalModel() const
{
    return gm_;
}

template<class GM, class ACCUMULATOR>
inline const size_t& 
LazyFlipper<GM, ACCUMULATOR>::maxSubgraphSize() const
{
    return maxSubgraphSize_;
}

template<class GM, class ACCUMULATOR>
inline void
LazyFlipper<GM, ACCUMULATOR>::setMaxSubgraphSize
(
    const size_t& maxSubgraphSize
)
{
    if(maxSubgraphSize < 1) {
        throw std::runtime_error("maximum subgraph size < 1.");
    }
    else {
        maxSubgraphSize_ = maxSubgraphSize;
    }
}

template<class GM, class ACCUMULATOR>
template<class VisitorType>
InferenceTermination
LazyFlipper<GM, ACCUMULATOR>::infer
(
    VisitorType& visitor
)
{
    size_t length = 1;
    for(;;) {
        visitor(*this, movemaker_.stateBegin(), movemaker_.stateEnd(), 
            movemaker_.energy(), length, subgraphForest_.size());
        SubgraphForestNode p = generateFirstPathOfLength(length);
        if(p == NONODE) {
            break;
        }
        else {
            while(p != NONODE) {
                if(Accumulator::bop(energyAfterFlip(p), movemaker_.energy())) {
                // if(energyAfterFlip(p) < movemaker_.energy()) {
                    flip(p);
                    activateInfluencedVariables(p, 0);
                    visitor(*this, movemaker_.stateBegin(), movemaker_.stateEnd(), movemaker_.energy(), 
                        length, subgraphForest_.size());
                }
                p = generateNextPathOfSameLength(p);
            }
            size_t currentActivationList = 0;
            size_t nextActivationList = 1;
            for(;;) {
                SubgraphForestNode p2 = firstActivePath(currentActivationList);
                if(p2 == NONODE) {
                    break;
                }
                else {
                    while(p2 != NONODE) {
                        if(Accumulator::bop(energyAfterFlip(p2), movemaker_.energy())) {
                        // if(energyAfterFlip(p2) < movemaker_.energy()) {
                            flip(p2);
                            activateInfluencedVariables(p2, nextActivationList);
                            visitor(*this, movemaker_.stateBegin(), movemaker_.stateEnd(), movemaker_.energy(), 
                                length, subgraphForest_.size());
                        }
                        p2 = nextActivePath(p2, currentActivationList);
                    }
                    deactivateAllVariables(currentActivationList);
                    nextActivationList = 1 - nextActivationList;
                    currentActivationList = 1 - currentActivationList;
                }
            }
        }
        if(length == maxSubgraphSize_) {
            break;
        }
        else {
            ++length;
        }
    }

    // assertion testing
    if(!NO_DEBUG) {
        subgraphForest_.testInvariant();
    }

    // diagnose
    // std::cout << subgraphForest_.asString();

    return NORMAL;
}

template<class GM, class ACCUMULATOR>
inline InferenceTermination
LazyFlipper<GM, ACCUMULATOR>::infer()
{
    LazyFlipperVisitor<LazyFlipper<gm_type, Accumulator> > v;
    return infer(v);
}

template<class GM, class ACCUMULATOR>
inline InferenceTermination 
LazyFlipper<GM, ACCUMULATOR>::arg
(
    std::vector<state_type>& arg, 
    const size_t& n
)const
{
    if(n > 1) {
        return UNKNOWN;
    }
    else {
        arg.resize(gm_.space().dimension());
        for(size_t j=0; j<gm_.space().dimension(); ++j) {
            arg[j] = movemaker_.state(j);
        }
        return NORMAL;
    }
}

template<class GM, class ACCUMULATOR>
inline typename LazyFlipper<GM, ACCUMULATOR>::value_type 
LazyFlipper<GM, ACCUMULATOR>::optimum() const
{
    return movemaker_.energy();
}

// Append the next possible variable to a node in the subgraph tree.
//
// The null pointer is returned if no variable can be appended.
// 
template<class GM, class ACCUMULATOR>
typename LazyFlipper<GM, ACCUMULATOR>::SubgraphForestNode
LazyFlipper<GM, ACCUMULATOR>::appendVariableToPath
(
    typename LazyFlipper<GM, ACCUMULATOR>::SubgraphForestNode p // input
)
{
    // collect variable indices on path
    std::vector<size_t> variableIndicesOnPath(subgraphForest_.level(p) + 1);
    {
        SubgraphForestNode p2 = p;
        for(size_t j=0; j<=subgraphForest_.level(p); ++j) {
            Assert(NO_DEBUG || p2 != NONODE);
            variableIndicesOnPath[subgraphForest_.level(p) - j] = subgraphForest_.value(p2);
            p2 = subgraphForest_.parent(p2);
        }
        Assert(NO_DEBUG || p2 == NONODE);
    }

    // find the mininum and maximum variable index on the path
    size_t minVI = variableIndicesOnPath[0];
    size_t maxVI = variableIndicesOnPath[0];
    for(size_t j=1; j<variableIndicesOnPath.size(); ++j) {
        if(variableIndicesOnPath[j] > maxVI) {
            maxVI = variableIndicesOnPath[j];
        }
    }
    // find the maximum variable index among the children of p.
    // the to be appended variable must have a greater index.
    if(subgraphForest_.numberOfChildren(p) > 0) {
        size_t maxChildIndex = subgraphForest_.numberOfChildren(p) - 1;
        minVI = subgraphForest_.value(subgraphForest_.child(p, maxChildIndex));
    }

    // build set of candidate variable indices for appending
    std::set<size_t> candidateVariableIndices;
    {
        SubgraphForestNode q = p;
        while(q != NONODE) {
            for(Adjacency::const_iterator it = variableAdjacency_.neighborsBegin(subgraphForest_.value(q));
                it != variableAdjacency_.neighborsEnd(subgraphForest_.value(q)); ++it) {
                    candidateVariableIndices.insert(*it);
            }
            q = subgraphForest_.parent(q);
        }
    }

    // append candidate if possible
    for(std::set<size_t>::const_iterator it = candidateVariableIndices.begin();
        it != candidateVariableIndices.end(); ++it) {
            // for all variables adjacenct to the one at node p
            if(*it > minVI && std::find(variableIndicesOnPath.begin(), variableIndicesOnPath.end(), *it) == variableIndicesOnPath.end()) {
                // the variable index *it is not smaller than the lower bound AND 
                // greater than the minimum variable index on the path AND 
                // is not itself on the path (??? consider tagging instead of 
                // searching in the previous if-condition)
                if(*it > maxVI) {
                    // *it is greater than the largest variable index on the path
                    return subgraphForest_.push_back(*it, p); // append to path
                }
                else { 
                    // *it is not the greatest variable index on the path.
                    for(size_t j=1; j<variableIndicesOnPath.size(); ++j) {
                        if(variableAdjacency_.connected(variableIndicesOnPath[j-1], *it)) {
                            // *it could have been added as a child of 
                            // variableIndicesOnPath[j-1]
                            for(size_t k=j; k<variableIndicesOnPath.size(); ++k) {
                                if(*it < variableIndicesOnPath[k]) {
                                    // adding *it as a child of variableIndicesOnPath[j-1]
                                    // would have made the path cheaper
                                    goto doNotAppend; // escape loop over j
                                }
                            }
                        }
                    }
                    // *it could not have been introduced cheaper
                    // append to path:
                    return subgraphForest_.push_back(*it, p); 
doNotAppend:;
                }
            }
    }
    // no neighbor of p could be appended
    return NONODE;
}

template<class GM, class ACCUMULATOR>
typename LazyFlipper<GM, ACCUMULATOR>::SubgraphForestNode 
LazyFlipper<GM, ACCUMULATOR>::generateFirstPathOfLength
(
    const size_t& length
)
{
    Assert(NO_ARG_TEST || length > 0);
    if(length > gm_.space().dimension()) {
        return NONODE;
    }
    else {
        if(length == 1) {
            SubgraphForestNode p = subgraphForest_.push_back(0, NONODE);
            // variable index = 0, parent = NONODE
            return p;
        }
        else {
            SubgraphForestNode p = subgraphForest_.levelAnchor(length-2);
            while(p != NONODE) {
                SubgraphForestNode p2 = appendVariableToPath(p);
                if(p2 != NONODE) { // append succeeded
                    return p2;
                }
                else { // append failed
                    p = subgraphForest_.levelOrderSuccessor(p);
                }
            }
            return NONODE;
        }
    }
}

template<class GM, class ACCUMULATOR>
typename LazyFlipper<GM, ACCUMULATOR>::SubgraphForestNode 
LazyFlipper<GM, ACCUMULATOR>::generateNextPathOfSameLength
(
    SubgraphForestNode predecessor
)
{
    if(subgraphForest_.level(predecessor) == 0) {
        if(subgraphForest_.value(predecessor) + 1 < gm_.space().dimension()) {
            SubgraphForestNode newNode = 
                subgraphForest_.push_back(subgraphForest_.value(predecessor) + 1, NONODE);
            subgraphForest_.setLevelOrderSuccessor(predecessor, newNode);
            return newNode;
        }
        else {
            // no more variables
            return NONODE;
        }
    }
    else {
        for(SubgraphForestNode parent = subgraphForest_.parent(predecessor);
            parent != NONODE; parent = subgraphForest_.levelOrderSuccessor(parent) ) {
                SubgraphForestNode newNode = appendVariableToPath(parent);
                if(newNode != NONODE) { 
                    // a variable has been appended
                    subgraphForest_.setLevelOrderSuccessor(predecessor, newNode);
                    return newNode;
                }
        }
        return NONODE;
    }
}

template<class GM, class ACCUMULATOR>
void 
LazyFlipper<GM, ACCUMULATOR>::activateInfluencedVariables
(
    SubgraphForestNode p,
    const size_t& activationListIndex
)
{
    Assert(NO_ARG_TEST || activationListIndex < 2);
    while(p != NONODE) {
        activation_[activationListIndex].tag(subgraphForest_.value(p), true);
        for(Adjacency::const_iterator it = variableAdjacency_.neighborsBegin(subgraphForest_.value(p));
            it != variableAdjacency_.neighborsEnd(subgraphForest_.value(p)); ++it) {
                activation_[activationListIndex].tag(*it, true);
        }
        p = subgraphForest_.parent(p);
    }
}

template<class GM, class ACCUMULATOR>
inline void 
LazyFlipper<GM, ACCUMULATOR>::deactivateAllVariables
(
    const size_t& activationListIndex
)
{
    Assert(NO_ARG_TEST || activationListIndex < 2);
    activation_[activationListIndex].untag();
}

template<class GM, class ACCUMULATOR>
typename LazyFlipper<GM, ACCUMULATOR>::SubgraphForestNode 
LazyFlipper<GM, ACCUMULATOR>::firstActivePath
(
    const size_t& activationListIndex
)
{
    if(subgraphForest_.levels() == 0) {
        return NONODE;
    }
    else {
        // ??? improve code: no search, store reference
        SubgraphForestNode p = subgraphForest_.levelAnchor(0);
        while(p != NONODE) {
            if(activation_[activationListIndex].tag(subgraphForest_.value(p))) {
                return p;
            }
            p = subgraphForest_.levelOrderSuccessor(p);
        }
        return NONODE;
    }
}

// ??? improve code: searching over all paths and all variables of each
// path for active variables is certainly not the ideal way
template<class GM, class ACCUMULATOR>
typename LazyFlipper<GM, ACCUMULATOR>::SubgraphForestNode 
LazyFlipper<GM, ACCUMULATOR>::nextActivePath
(
    SubgraphForestNode predecessor, 
    const size_t& activationListIndex
)
{
    for(;;) {
        if(subgraphForest_.levelOrderSuccessor(predecessor) == NONODE) {
            if(subgraphForest_.level(predecessor) + 1 < subgraphForest_.levels()) {
                // there are more levels in the tree
                predecessor = subgraphForest_.levelAnchor(subgraphForest_.level(predecessor) + 1);
            }
            else {
                // there are no more levels in the tree
                return NONODE;
            }
        }
        else {
            // predecessor is not the last node on its level
            predecessor = subgraphForest_.levelOrderSuccessor(predecessor);
        }
        SubgraphForestNode p = predecessor;
        while(p != NONODE) {
            // search along path for active variables:
            if(activation_[activationListIndex].tag(subgraphForest_.value(p))) {
                return predecessor;
            }
            p = subgraphForest_.parent(p);
        }
    }
}

template<class GM, class ACCUMULATOR>
inline typename LazyFlipper<GM, ACCUMULATOR>::value_type 
LazyFlipper<GM, ACCUMULATOR>::energyAfterFlip
(
    SubgraphForestNode node
)
{
    size_t numberOfFlippedVariables = subgraphForest_.level(node) + 1;
    std::vector<size_t> flippedVariableIndices(numberOfFlippedVariables);
    std::vector<state_type> flippedVariableStates(numberOfFlippedVariables);
    for(size_t j=0; j<numberOfFlippedVariables; ++j) {
        Assert(NO_DEBUG || node != NONODE);
        flippedVariableIndices[j] = subgraphForest_.value(node);
        // binary flip:
        flippedVariableStates[j] = 1 - movemaker_.state(subgraphForest_.value(node)); 
        node = subgraphForest_.parent(node);
    }
    Assert(NO_DEBUG || node == NONODE);
    return movemaker_.energyAfterMove(flippedVariableIndices.begin(), 
        flippedVariableIndices.end(), flippedVariableStates.begin());
}

template<class GM, class ACCUMULATOR>
inline void
LazyFlipper<GM, ACCUMULATOR>::flip
(
    SubgraphForestNode node
)
{
    size_t numberOfFlippedVariables = subgraphForest_.level(node) + 1;
    std::vector<size_t> flippedVariableIndices(numberOfFlippedVariables);
    std::vector<state_type> flippedVariableStates(numberOfFlippedVariables);
    for(size_t j=0; j<numberOfFlippedVariables; ++j) {
        Assert(NO_DEBUG || node != NONODE);
        flippedVariableIndices[j] = subgraphForest_.value(node);
        // binary flip:
        flippedVariableStates[j] = 1 - movemaker_.state(subgraphForest_.value(node)); 
        node = subgraphForest_.parent(node);
    }
    Assert(NO_DEBUG || node == NONODE);
    movemaker_.move(flippedVariableIndices.begin(), 
        flippedVariableIndices.end(), flippedVariableStates.begin());
}

// implementatin of Visitors

template<class LAZYFLIPPER>
template<class StateIterator>
inline void
LazyFlipperVisitor<LAZYFLIPPER>::operator()
(
    const typename LazyFlipperVisitor<LAZYFLIPPER>::lazy_flipper_type& lazyFlipper,
    StateIterator begin,
    StateIterator end,
    const typename LazyFlipperVisitor<LAZYFLIPPER>::value_type& energy,
    const size_t& subgraphSize,
    const size_t& subgraphForestSize
) const
{}

template<class LAZYFLIPPER>
inline 
LazyFlipperVerboseVisitor<LAZYFLIPPER>::LazyFlipperVerboseVisitor()
: step_(0)
{}

template<class LAZYFLIPPER>
template<class StateIterator>
inline void
LazyFlipperVerboseVisitor<LAZYFLIPPER>::operator()
(
    const typename LazyFlipperVerboseVisitor<LAZYFLIPPER>::lazy_flipper_type& lazyFlipper,
    StateIterator begin,
    StateIterator end,
    const typename LazyFlipperVerboseVisitor<LAZYFLIPPER>::value_type& energy,
    const size_t& subgraphSize,
    const size_t& subgraphForestSize
)
{
    ++step_;
    std::cout << "step " << step_ 
        << ": E=" << energy
        << ", subgraph_size=" << subgraphSize
        << ", subgraph_tree_size=" << subgraphForestSize
        << std::endl;
}

} // namespace opengm

#endif // #ifndef OPENGM_LAZYFLIPPER_HXX
