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
#ifndef OPENGM_DECOMPOSER_HXX
#define OPENGM_DECOMPOSER_HXX

#include <vector>
#include <list>
#include <map>

#include <ufd/ufd.hxx>

#include "opengm/graphicalmodel.hxx"

namespace opengm {

template<class GM>
class Decomposer {
public:
    typedef GM GraphicalModel;
    typedef typename GraphicalModel::factor_type Factor;

    void decompose(const GraphicalModel&, std::vector<std::vector<size_t> >&);
};

template<class GM>
void
Decomposer<GM>::decompose
(
    const GraphicalModel& gm,
    std::vector<std::vector<size_t> >& decomposition
)
{
    std::list<size_t> blackList;
    std::list<size_t> grayList;
    std::list<size_t> whiteList;
    for(size_t j=0; j<gm.numberOfFactors(); ++j) {
        whiteList.push_back(j);
    }

    size_t n = 0;
    while(!whiteList.empty()) {
        decomposition.resize(n+1);
        ufd::Partition<size_t> partition(gm.space().dimension());
        std::list<size_t>* currentList;
        for(size_t listSwap=0; listSwap<2; ++listSwap) {
            if(listSwap == 0) {
                currentList = &whiteList;
            }
            else {
                currentList = &blackList;
            }

            std::list<size_t>::iterator it = currentList->begin();
            while(it != currentList->end()) {
                // check if *it can be inserted in decomposition[n]
                bool insert = true;
                const Factor& factor = gm[*it];
                std::map<size_t, size_t> counters;
                for(size_t j=0; j<factor.numberOfVariables(); ++j) {
                    size_t c = ++counters[partition.find(factor.variableIndex(j))];
                    if(c > 1) {
                        insert = false;
                        break;
                    }
                }
                if(insert) {
                    decomposition[n].push_back(*it);
                    if(currentList == &whiteList) {
                        grayList.push_back(*it);
                        it = currentList->erase(it);
                    }
                    else {
                        ++it;
                    }
                    for(size_t j=1; j<factor.numberOfVariables(); ++j) {
                        partition.merge(factor.variableIndex(j-1), factor.variableIndex(j));
                    }
                }
                else {
                    ++it;
                }
            }
        }
        blackList.insert(blackList.end(), grayList.begin(), grayList.end());
        grayList.clear();
        ++n;
    }
}

} // namespace opengm

#endif // #idndef OPENGM_DECOMPOSER_HXX
