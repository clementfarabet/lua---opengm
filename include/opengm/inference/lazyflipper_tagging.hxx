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
#ifndef OPENGM_LAZYFLIPPER_TAGGING_HXX
#define OPENGM_LAZYFLIPPER_TAGGING_HXX

#include <vector>
#include <list>

namespace opengm {

template<class T>
class Tagging
{
public:
    typedef T value_type;
    typedef std::vector<value_type> tag_container_type;
    typedef std::vector<size_t> index_container_type;
    typedef index_container_type::const_iterator const_iterator;

    Tagging(const size_t& = 0);

    void append(const size_t&);
    value_type tag(const size_t&) const;
    void tag(const size_t&, const value_type&);
    void untag(); // untag all
    const_iterator begin();
    const_iterator begin() const;
    const_iterator end();
    const_iterator end() const;

private:
    tag_container_type tags_;
    index_container_type indices_;
};

// implementation of Tagging

template<class T>
inline Tagging<T>::Tagging
(
    const size_t& size
)
:   tags_(tag_container_type(size)), 
    indices_(index_container_type())
{}

template<class T>
inline void Tagging<T>::append
(
    const size_t& number
)
{
    tags_.resize(tags_.size() + number);
}

// runtime complexity: constant
template<class T>
inline typename Tagging<T>::value_type
Tagging<T>::tag
(
    const size_t& index
) const
{
    Assert(NO_DEBUG || index < tags_.size());
    return tags_[index];
}

// runtime complexity: constant
template<class T>
inline void 
Tagging<T>::tag
(
    const size_t& index,
    const typename Tagging<T>::value_type& tag
)
{
    Assert(NO_DEBUG || index < tags_.size());
    Assert(NO_DEBUG || tag != T()); // no implicit un-tagging
    if(tags_[index] == T()) { // so far un-tagged
        indices_.push_back(index);
    }
    tags_[index] = tag;
}

// untag all
// runtime complexity: linear in indices_.size()
// note the performance gain over linearity in tags_.size()
template<class T>
inline void 
Tagging<T>::untag()
{
    for(const_iterator it = indices_.begin(); it != indices_.end(); ++it) {
        tags_[*it] = T();
    }
    indices_.clear();
}

template<class T>
inline typename Tagging<T>::const_iterator
Tagging<T>::begin() const
{
    return indices_.begin();
}

template<class T>
inline typename Tagging<T>::const_iterator
Tagging<T>::end() const
{
    return indices_.end();
}

template<class T>
inline typename Tagging<T>::const_iterator
Tagging<T>::begin() 
{
    return indices_.begin();
}

template<class T>
inline typename Tagging<T>::const_iterator
Tagging<T>::end()
{
    return indices_.end();
}

} // namespace opengm

#endif // #ifndef OPENGM_LAZYFLIPPER_TAGGING_HXX
