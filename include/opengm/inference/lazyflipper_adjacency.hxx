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
#ifndef OPENGM_LAZYFLIPPER_ADJACENCY_HXX
#define OPENGM_LAZYFLIPPER_ADJACENCY_HXX

namespace opengm {

// A simple undirected graph
class Adjacency {
public:
    typedef std::set<size_t>::const_iterator const_iterator;

    Adjacency(const size_t& = 0);
    void resize(const size_t&);
    void connect(const size_t&, const size_t&);
    bool connected(const size_t&, const size_t&) const;
    const_iterator neighborsBegin(const size_t&);
    const_iterator neighborsBegin(const size_t&) const;
    const_iterator neighborsEnd(const size_t&);
    const_iterator neighborsEnd(const size_t&) const;

private:
    std::vector<std::set<size_t> > neighbors_;
};

// implementation of Adjacency

inline
Adjacency::Adjacency
(
    const size_t& size
) : neighbors_(std::vector<std::set<size_t> >(size))
{}

inline void 
Adjacency::resize
(
    const size_t& size
)
{
    neighbors_.resize(size);
}

inline void 
Adjacency::connect
(
    const size_t& j, 
    const size_t& k
)
{
    neighbors_[j].insert(k);
    neighbors_[k].insert(j);
}

inline bool 
Adjacency::connected
(
    const size_t& j, 
    const size_t& k
) const
{
    if(neighbors_[j].size() < neighbors_[k].size()) {
        if(neighbors_[j].find(k) == neighbors_[j].end()) {
            return false;
        }
        else {
            return true;
        }
    }
    else {
        if(neighbors_[k].find(j) == neighbors_[k].end()) {
            return false;
        }
        else {
            return true;
        }
    }
}

inline Adjacency::const_iterator 
Adjacency::neighborsBegin
(
    const size_t& index
)
{
    return neighbors_[index].begin();
}

inline Adjacency::const_iterator 
Adjacency::neighborsBegin
(
    const size_t& index
) const
{
    return neighbors_[index].begin();
}

inline Adjacency::const_iterator 
Adjacency::neighborsEnd
(
    const size_t& index
)
{
    return neighbors_[index].end();
}

inline Adjacency::const_iterator 
Adjacency::neighborsEnd
(
    const size_t& index
) const
{
    return neighbors_[index].end();
}

} // namespace opengm

#endif // #ifndef OPENGM_LAZYFLIPPER_ADJACENCY_HXX
