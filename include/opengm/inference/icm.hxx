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
#ifndef OPENGM_ICM_HXX
#define OPENGM_ICM_HXX

#include <vector>
#include <string>
#include <iostream>

#include "opengm/opengm.hxx"
#include "opengm/graphicalmodel.hxx"
#include "opengm/inference/inference.hxx"
#include "opengm/inference/movemaker.hxx"

namespace opengm {

template<class GM, class ACCUMULATOR>
class ICM : public Inference<GM, ACCUMULATOR>
{
public:	    
    typedef GM                                  gm_type;
    typedef ACCUMULATOR                         Accumulator;
    typedef typename gm_type::factor_type       factor_type;
    typedef typename factor_type::value_type    value_type;   
    typedef typename factor_type::space_type    space_type;
    typedef typename space_type::state_type     state_type;
    typedef typename gm_type::Operator          Operator;
    typedef Movemaker<gm_type>                  movemaker_type;

    class Parameter{
    public:
        std::vector<size_t> startPoint_;
    };

    ICM(const gm_type&);
    ICM(const gm_type&, const Parameter&);
    virtual std::string name() const;
    virtual const gm_type& graphicalModel() const;
    virtual InferenceTermination infer();
    template<class VisitorType>
        InferenceTermination infer(VisitorType&);
    virtual InferenceTermination arg(std::vector<state_type>&, const size_t& = 1)const ;

private:
    const gm_type& gm_;  
    movemaker_type movemaker_;
}; 

// visitors

template<class ICM>
class ICMVisitor {
public:
    typedef ICM icm_type;
    typedef typename icm_type::value_type value_type;

    template<class StateIterator>
        void operator()(const icm_type&, StateIterator, StateIterator, const value_type&) const;
};

template<class ICM>
class ICMVerboseVisitor {
public:
    typedef ICM icm_type;
    typedef typename icm_type::value_type value_type;

    ICMVerboseVisitor();
    template<class StateIterator>
        void operator()(const icm_type&, StateIterator, StateIterator, const value_type&);

private:
    size_t step_;
};

// implementation of ICM

template<class GM, class ACCUMULATOR>
ICM<GM, ACCUMULATOR>::ICM
(
    const gm_type& gm
) : gm_(gm), movemaker_(gm)
{
} 

template<class GM, class ACCUMULATOR>
ICM<GM, ACCUMULATOR>::ICM
(
    const gm_type& gm,
    const Parameter& parameter
) : gm_(gm), movemaker_(gm, parameter.startPoint_.begin())
{
} 

template<class GM, class ACCUMULATOR>
inline std::string 
ICM<GM, ACCUMULATOR>::name() const
{
    return "ICM";
}

template<class GM, class ACCUMULATOR>
inline const typename ICM<GM, ACCUMULATOR>::gm_type&
ICM<GM, ACCUMULATOR>::graphicalModel() const
{
    return gm_;
} 

template <class GM, class ACCUMULATOR>
inline InferenceTermination 
ICM<GM,ACCUMULATOR>::infer()
{
    ICMVisitor<ICM<gm_type, Accumulator> > v;
    return infer(v);
}

template<class GM, class ACCUMULATOR>
template<class VisitorType>
InferenceTermination ICM<GM,ACCUMULATOR>::infer
(
    VisitorType& visitor
)
{
    visitor(*this, movemaker_.stateBegin(), movemaker_.stateEnd(), movemaker_.energy());
    bool updates = true;
    while(updates) {
        updates = false;
        for(size_t v=0; v<gm_.space().dimension(); ++v) {
            for(size_t s=0; s<gm_.space().numberOfStates(v); ++s) {
                if(s != movemaker_.state(v)) {
                    if(Accumulator::bop(movemaker_.energyAfterMove(&v, &v+1, &s), movemaker_.energy())) {
                        movemaker_.move(&v, &v+1, &s);
                        updates = true;
                        visitor(*this, movemaker_.stateBegin(), movemaker_.stateEnd(), movemaker_.energy());
                    }
                }
            }
        }
    }
    return CONVERGENCE;
}

template <class GM, class ACCUMULATOR>
inline InferenceTermination 
ICM<GM,ACCUMULATOR>::arg
(
    std::vector<state_type>& x, 
    const size_t& N
)const
{
    if(N==1){
        x.resize(gm_.space().dimension());
        for(size_t j=0; j<x.size(); ++j) {
            x[j] = movemaker_.state(j);
        }
        return NORMAL;
    }
    else{
        return UNKNOWN;
    }
}

// implementation of visitors

template<class ICM>
template<class StateIterator>
inline void
ICMVisitor<ICM>::operator()
(
    const typename ICMVisitor<ICM>::icm_type& icm,
    StateIterator stateBegin,
    StateIterator stateEnd,
    const typename ICMVisitor<ICM>::value_type& value
) const
{}

template<class ICM>
inline 
ICMVerboseVisitor<ICM>::ICMVerboseVisitor()
: step_(0)
{
}

template<class ICM>
template<class StateIterator>
inline void
ICMVerboseVisitor<ICM>::operator()
(
    const typename ICMVerboseVisitor<ICM>::icm_type& icm,
    StateIterator stateBegin,
    StateIterator stateEnd,
    const typename ICMVerboseVisitor<ICM>::value_type& value
)
{
    ++step_;
    std::cout << "step " << step_ << ": E=" << value << std::endl;
}

} // namespace opengm

#endif
