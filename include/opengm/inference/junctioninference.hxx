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
#ifndef OPENGM_JUNCTIONINFERNCE_HXX
#define OPENGM_JUNCTIONINFERNCE_HXX

#include <vector>
#include <iostream>

#include <marray/marray.hxx>

#include "opengm/opengm.hxx"
#include "opengm/inference/inference.hxx"

namespace opengm
{

template <class GM, class ACC, class INF>
class JunctionInference : public Inference<GM,ACC>
{
public:
    typedef typename GM::factor_type          factor_type;
    typedef typename factor_type::value_type  value_type;   
    typedef typename factor_type::space_type  space_type;
    typedef typename space_type::state_type   state_type;
    typedef typename GM::Operator             OP;
    typedef typename std::set<state_type>     container;

    class Parameter{
    public:
        std::vector<container > jNodeIndices_;
        std::vector<container > jFactorJIndices_;
        typename INF::Parameter parameter_;
    };   

    // construction and destruction
    JunctionInference(const GM& gm, Parameter para=Parameter());
    ~JunctionInference();

    // inference
    virtual InferenceTermination infer();
    virtual InferenceTermination arg(std::vector<state_type>&, const size_t& = 1);

    // query
    virtual std::string name() const 
    {return "JunctionInference";} 
    const GM& graphicalModel() const;
    const GM& graphicalJunctionModel() const;

private:
    void mapOSpace2JSpace(const std::vector<size_t>& oconf,
        const std::vector<size_t>& jvarsize,
        const std::vector<size_t>& var,  
        const std::vector<size_t>& jmap, 
        std::vector<size_t>& jconf) const;
    void mapJSpace2OSpace(const std::vector<size_t>& jconf, 
        const std::vector<size_t>& jvar, 
        std::vector<size_t>& oconf) const;

    const GM& gm_;
    GM jgm_;
    DiscreteSpace jspace_;
    Parameter para_;
    INF* inference_; 
    //std::set<std::set<state_type> >    jFactorIndices_;
    //std::vector<std::vector<state_type> > jFactorIndicesJ_;
};


template <class GM, class ACC, class INF>
JunctionInference<GM,ACC,INF>::JunctionInference
(
    const GM& gm,
    Parameter para
) : gm_(gm)
{
    para_= para;
    size_t numJNodes   = para_.jNodeIndices_.size();
    size_t numJFactors = para_.jFactorJIndices_.size();
    //Compute number of states of each junction-graph-variable 
    std::vector<size_t> numberOfStates(numJNodes,1);
    for(size_t jnode=0; jnode<numJNodes; ++jnode){
        typename container::iterator it;
        for(it=para_.jNodeIndices_[jnode].begin(); it!=para_.jNodeIndices_[jnode].end(); ++it){
            numberOfStates[jnode] *= gm_.space().numberOfStates(*it);
        }
    } 

    //Build Junction-Space & Junction-GM 
    jspace_ = DiscreteSpace(numberOfStates.begin(), numberOfStates.end());
    for(size_t jfac=0; jfac<numJFactors; ++jfac){
        //allocate jfactor in jspace
        assert(para_.jFactorJIndices_[jfac].size()>0);
        typename container::iterator it;
        for(it=para_.jFactorJIndices_[jfac].begin(); it !=para_.jFactorJIndices_[jfac].end(); ++it){
            assert(*it<numJNodes);
        }
        jgm_.addFactor(factor_type(jspace_,
            para_.jFactorJIndices_[jfac].begin(),
            para_.jFactorJIndices_[jfac].end(),
            ACC::template neutral<value_type>()
            )
            );
    }

    //Copy data into the junction graph 
    std::vector<bool> factorHasBeenAdded(gm_.numberOfFactors(),false);
    for(size_t jfac=0; jfac<numJFactors; ++jfac){
        //Build factor in original space
        std::set<size_t> var;
        std::vector<size_t> jvar;
        typename container::iterator it; 
        typename container::iterator it2;
        for(it=para_.jFactorJIndices_[jfac].begin() ; it!=para_.jFactorJIndices_[jfac].end(); ++it){
            for(it2=para_.jNodeIndices_[*it].begin() ; it2!=para_.jNodeIndices_[*it].end() ; ++it2){
                var.insert(*it2);
                jvar.push_back(*it2);
            }
        }
        factor_type factor(gm_.space(), var.begin(), var.end(),OP::template neutral<value_type>());

        for(size_t fac=0; fac<gm_.numberOfFactors(); ++fac){
            if( factorHasBeenAdded[fac] ) continue;
            bool covers = true;
            for(size_t i=0; i<gm_[fac].numberOfVariables(); ++i){
                if(var.count(gm_[fac].variableIndex(i))==0){
                    covers = false;
                    continue; // jfac does not cover fac
                }
            }
            if(covers){
                //invariant: jfac coves fac and fac has not been added
                OP::op(gm_[fac], factor);
                factorHasBeenAdded[fac] = true;
            }
        }

        //Copy data from factor to jgm[jfac]
        std::vector<size_t> seq(factor.table().dimension());
        std::vector<size_t> jvarsize(para_.jFactorJIndices_[jfac].size());
        std::vector<size_t> vecvar(var.size());
        std::vector<size_t> jconf(jgm_[jfac].numberOfVariables());
        std::vector<size_t> jmap(jvar.size());
        size_t count;
        for(size_t i=0; i<jmap.size(); ++i){
            count = 0;
            for(std::set<size_t>::iterator it= var.begin(); it!=var.end();++it){
                if(*it==jvar[i]){
                    jmap[i] = count;
                    break;
                }
                ++count;
            }
        }
        count=0;
        for(std::set<size_t>::iterator it= var.begin(); it!=var.end();++it){
            vecvar[count++]=*it;
        }
        count=0;
        for(std::set<size_t>::iterator it=para_.jFactorJIndices_[jfac].begin(); it!=para_.jFactorJIndices_[jfac].end();++it){
            jvarsize[count++]=para_.jNodeIndices_[*it].size();
        }
        for(typename marray::Marray<value_type>::const_iterator it = factor.table().begin(); it.hasMore(); ++it) {
            it.coordinate(seq.begin());
            mapOSpace2JSpace(seq,jvarsize,vecvar,jmap,jconf);
            jgm_[jfac](jconf.begin()) = *it;
        }

    } 
    for(size_t fac=0; fac<gm_.numberOfFactors(); ++fac){
        assert(factorHasBeenAdded[fac]==true);
    }

    inference_ = new INF(jgm_,para_.parameter_);
}
template <class GM, class ACC, class INF>
JunctionInference<GM,ACC,INF>::~JunctionInference()
{
    delete inference_;
}

template <class GM, class ACC, class INF>
void JunctionInference<GM,ACC,INF>::mapOSpace2JSpace
(
    const std::vector<size_t>& oconf, 
    const std::vector<size_t>& jvarsize,
    const std::vector<size_t>& var,
    const std::vector<size_t>& jmap,
    std::vector<size_t>& jconf
)const
{
    size_t c = 0;
    for(size_t i=0; i<jvarsize.size(); ++i){
        jconf[i] = oconf[jmap[c++]];
        for(size_t j=1; j<jvarsize[i]; ++j){
            jconf[i] = jconf[i] * gm_.space().numberOfStates(var[c]) + oconf[jmap[c]];
            c++;
        }
    }
}

template <class GM, class ACC, class INF>
void JunctionInference<GM,ACC,INF>::mapJSpace2OSpace
(
    const std::vector<size_t>& jconf,
    const std::vector<size_t>& jvar, 
    std::vector<size_t>& oconf
) const
{
    assert(oconf.size()==gm_.space().dimension());
    for(size_t i=0; i<jvar.size(); ++i){
        std::vector<size_t> numStates(para_.jNodeIndices_[i].size(),1);
        std::set<size_t>::const_iterator it;

        it = para_.jNodeIndices_[i].begin();
        for(size_t j=0; j<para_.jNodeIndices_[i].size(); ++j){
            numStates[j] = gm_.space().numberOfStates(*it);	
            ++it;
        }

        //size_t aa=5; //FIXME unused
        //size_t bb=2; //FIXME unsued
        assert(5/2==2);
        assert(5%2==1);



        size_t ind = jconf[i];
        std::vector<size_t> t(para_.jNodeIndices_[i].size());
        size_t n = para_.jNodeIndices_[i].size();
        do{
            --n;
            t[n] = ind%numStates[n];
            ind = ind/numStates[n];
        }while(n!=0);

        it = para_.jNodeIndices_[i].begin();
        for(size_t j=0; j<para_.jNodeIndices_[i].size(); ++j){
            oconf[*it] = t[j];
            ++it;
        }
    }
}

template <class GM, class ACC, class INF>
const GM& JunctionInference<GM,ACC,INF>::graphicalModel() const
{
    return gm_;
}

template <class GM, class ACC, class INF>
const GM& JunctionInference<GM,ACC,INF>::graphicalJunctionModel() const
{
    return jgm_;
}

template <class GM, class ACC, class INF>
InferenceTermination JunctionInference<GM,ACC,INF>::infer()
{
    return inference_->infer();
}

template <class GM, class ACC, class INF>
InferenceTermination JunctionInference<GM,ACC,INF>::arg
(
    std::vector<state_type>& conf, 
    const size_t& n
)
{
    conf.resize(gm_.space().dimension());
    std::vector<state_type> jconf;
    std::vector<size_t> jvar(jgm_.space().dimension());
    for(size_t i=0; i<jgm_.space().dimension(); ++i){
        jvar[i] = i; 
    }
    InferenceTermination r = inference_->arg(jconf,n);    
    mapJSpace2OSpace(jconf,jvar,conf);   
    return r;
}

} // namespace opengm

#endif // #ifndef OPENGM_JUNCTIONINFERNCE_HXX
