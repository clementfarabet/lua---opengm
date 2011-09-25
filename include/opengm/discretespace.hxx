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
#ifndef OPENGM_DISCRETESPACE_HXX
#define OPENGM_DISCRETESPACE_HXX

#include <vector>

#include "opengm.hxx"

namespace opengm {

class DiscreteSpace {
public: 
    typedef size_t state_type;
    typedef std::vector<size_t> container_type;

    DiscreteSpace();
    template<class Iterator>
        DiscreteSpace(Iterator, Iterator);
    template<class Iterator>
        void assign(Iterator, Iterator);

    size_t dimension() const;
    size_t numberOfStates(size_t) const;
    const container_type& numbersOfStates() const;
    template<class InIterator, class OutIterator>
        void numbersOfStates(InIterator, InIterator, OutIterator) const;
    template<class IndexIterator>
        bool isValidIndexSequence(IndexIterator, IndexIterator) const;

private:
    container_type numbersOfStates_;
};

inline
DiscreteSpace::DiscreteSpace()
: numbersOfStates_(container_type())
{}

template<class Iterator>
inline 
DiscreteSpace::DiscreteSpace
(
    Iterator begin,
    Iterator end
)
{
    assign(begin, end);
}

template<class Iterator>
inline void
DiscreteSpace::assign
(
    Iterator begin,
    Iterator end
)
{
    Assert(NO_ARG_TEST || begin != end);
    size_t dimension = std::distance(begin, end);
    container_type nos(dimension);
    container_type::iterator it = nos.begin();
    while(begin != end)	{
        Assert(NO_ARG_TEST || *begin != 0);
        *it = *begin;
        ++it;
        ++begin;
    }
    numbersOfStates_ = nos;
}

inline size_t
DiscreteSpace::dimension() const
{
    return numbersOfStates_.size();
}

inline size_t
DiscreteSpace::numberOfStates
(
    size_t dimension
) const
{
    Assert(NO_ARG_TEST || dimension < numbersOfStates_.size());
    return numbersOfStates_[dimension];
}

template<class InIterator, class OutIterator>
inline void
DiscreteSpace::numbersOfStates
(
    InIterator begin,
    InIterator end,
    OutIterator out
) const
{
    while(begin != end) {
        Assert(NO_ARG_TEST || *begin < numbersOfStates_.size());
        *out = numbersOfStates_[size_t(*begin)];
        ++begin;
        ++out;
    }
}

inline const DiscreteSpace::container_type& 
DiscreteSpace::numbersOfStates() const
{
    return numbersOfStates_;
}

// test a sequence of variable indices for monotonicity and bounds 
template<class IndexIterator>
inline bool
DiscreteSpace::isValidIndexSequence
(
    IndexIterator begin,
    IndexIterator end
) const
{
    IndexIterator previousIt = begin;
    while(begin != end) {
        if(*begin >= this->dimension()) {
            return false;
        }
        if(previousIt != begin && *previousIt >= *begin) {
            return false;
        }
        previousIt = begin;
        ++begin;
    }
    return true;
}

} // namespace opengm

#endif // #ifndef OPENGM_DISCRETESPACE_HXX
