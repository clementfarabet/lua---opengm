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
#ifndef OPENGM_SERIALIZATION_HXX
#define OPENGM_SERIALIZATION_HXX

#include "opengm/explicitfactor.hxx"
#include "opengm/graphicalmodel.hxx"

namespace opengm {

template<class GraphicalModel, class U>
void serialize
(
    const GraphicalModel& gm,
    marray::Vector<U>& out
)
{
    // compute size of serialization
    size_t size = 2; // number of variables, number of factors
    size += gm.space().dimension(); // number of states for each variable
    for(size_t j=0; j<gm.numberOfFactors(); ++j) {
        const typename GraphicalModel::factor_type& factor = gm[j];
        size += 1; // number of variables
        size += factor.numberOfVariables(); // variable indices
        size += factor.table().size(); // value table
    }
    out.resize(size);

    // serialize
    out[0] = static_cast<U>(gm.space().dimension());
    out[1] = static_cast<U>(gm.numberOfFactors());
    size_t p = 1;
    for(size_t j=0; j<gm.space().dimension(); ++j) {
        ++p;
        out[p] = static_cast<U>(gm.space().numberOfStates(j));
    }
    for(size_t j=0; j<gm.numberOfFactors(); ++j) {
        const typename GraphicalModel::factor_type& factor = gm[j];
        ++p;
        out[p] = static_cast<U>(factor.numberOfVariables());
        for(size_t k=0; k<factor.numberOfVariables(); ++k) {
            ++p;
            out[p] = static_cast<U>(factor.variableIndex(k));
        }
        for(size_t k=0; k<factor.table().size(); ++k) {
            ++p;
            out[p] = factor.table()(k);
        }
    }
    Assert(NO_DEBUG || p == out.size()-1);
}

template<class U, bool isConst, class GraphicalModel>
void deserialize
(
    const marray::View<U, isConst>& in,
    GraphicalModel& gm, 
    typename GraphicalModel::space_type& space
)
{
    size_t numberOfVariables = static_cast<size_t>(in(0));
    size_t numberOfFactors = static_cast<size_t>(in(1));
    std::vector<size_t> numbersOfStates(numberOfVariables);
    size_t p = 1;
    for(size_t j=0; j<numberOfVariables; ++j) {
        ++p;
        numbersOfStates[j] = static_cast<size_t>(in(p));
    }
    space = DiscreteSpace(numbersOfStates.begin(), numbersOfStates.end());

    gm = GraphicalModel(numberOfFactors);
    for(size_t j=0; j<numberOfFactors; ++j) {
        ++p;
        size_t numberOfVariablesOfFactor = static_cast<size_t>(in(p));
        std::vector<size_t> variableIndicesOfFactor(numberOfVariablesOfFactor);
        size_t tableSize = 1;
        for(size_t k=0; k<numberOfVariablesOfFactor; ++k) {
            ++p;
            variableIndicesOfFactor[k] = static_cast<size_t>(in(p));
            tableSize *= space.numberOfStates(variableIndicesOfFactor[k]);
        }
        std::vector<typename GraphicalModel::value_type> table(tableSize);
        for(size_t k=0; k<tableSize; ++k) {
            ++p;
            table[k] = in(p);
        }
        typename GraphicalModel::factor_type factor(space, variableIndicesOfFactor.begin(), 
            variableIndicesOfFactor.end(), table.begin(), table.end());
        gm.addFactor(factor);
    }
    Assert(NO_DEBUG || p == in.size()-1);
}

} // namespace opengm

#endif // #ifndef OPENGM_SERIALIZATION_HXX
