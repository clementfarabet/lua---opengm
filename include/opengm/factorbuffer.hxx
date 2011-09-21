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
#ifndef OPENGM_FACTORBUFFER_HXX
#define OPENGM_FACTORBUFFER_HXX

#include <vector>
#include <map>
#include <list>
#include <set>

namespace opengm {

template<class Factor>
class FactorBuffer
{ 
public: 
    typedef Factor factor_type;
    typedef typename Factor::space_type space_type;   
    typedef typename Factor::value_type value_type;

    // construction and assignment
    FactorBuffer();
    template<class IndexIterator>
        FactorBuffer(const space_type&, IndexIterator, IndexIterator,
            const value_type& = value_type());
    template<class IndexIterator>
        void assign(const space_type&, IndexIterator, IndexIterator,
            const value_type& = value_type());

    // access and manipulation
    Factor& current();
    Factor& old();
    const Factor& current() const;
    const Factor& old() const;
    void toggle();

    // distance between current and old
    template<class DIST> value_type dist() const;

private:
    bool flag_;
    Factor buffer1_;
    Factor buffer2_;
};

template<class Factor>
inline 
FactorBuffer<Factor>::FactorBuffer()
{}

template<class Factor>
template<class IndexIterator>
inline 
FactorBuffer<Factor>::FactorBuffer
(
    const space_type& space, 
    IndexIterator begin, 
    IndexIterator end,
    const typename Factor::value_type& constant
)
{
    assign(space, begin, end, constant);
} 

template<class Factor>
template<class IndexIterator>
inline void 
FactorBuffer<Factor>::assign
(
    const space_type& space, 
    IndexIterator begin, 
    IndexIterator end,
    const typename Factor::value_type& constant
)
{
    if(begin == end) {
        buffer1_.assign(space, constant);
        buffer2_.assign(space, constant);
    }
    else {
        buffer1_.assign(space, begin, end, constant);
        buffer2_.assign(space, begin, end, constant);
    } 
    flag_ = false;
}

template<class Factor>
inline Factor& 
FactorBuffer<Factor>::current() 
{
    return flag_ ? buffer1_ : buffer2_; 
};

template<class Factor>
inline const Factor& 
FactorBuffer<Factor>::current() const 
{
    return flag_ ? buffer1_:buffer2_; 
};

template<class Factor>
inline Factor& 
FactorBuffer<Factor>::old() 
{
    return flag_ ? buffer2_ : buffer1_; 
};

template<class Factor>
inline const Factor& 
FactorBuffer<Factor>::old() const 
{
    return flag_ ? buffer2_:buffer1_; 
};

template<class Factor>
inline void 
FactorBuffer<Factor>::toggle() 
{
    flag_ = flag_ ? false : true; 
};

template<class Factor>
template<class DIST> 
inline typename Factor::value_type 
FactorBuffer<Factor>::dist() const
{
    return DIST::op(buffer1_, buffer2_);
}

} // namespace opengm

#endif // #ifndef OPENGM_FACTORBUFFER_HXX
