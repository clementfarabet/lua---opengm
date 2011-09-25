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
#ifndef OPENGM_EXPLICITFACTOR_HXX
#define OPENGM_EXPLICITFACTOR_HXX

#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <marray/marray.hxx>

#include "opengm.hxx"
#include "discretespace.hxx"
#include "accumulation.hxx"
#include "minimizer.hxx"
#include "integrator.hxx"
#include "maximizer.hxx"

namespace opengm {

template<class T, class TABLE = marray::Marray<T> >
class ExplicitFactor {
public:
    typedef T value_type;
    typedef T& reference;
    typedef T* pointer;
    typedef const T& const_reference;
    typedef const T* const_pointer;
    typedef DiscreteSpace space_type;
    typedef typename space_type::state_type state_type;
    typedef unsigned int variable_index_type;
    typedef std::vector<variable_index_type> variable_index_container_type;
    typedef TABLE table_type;

    ExplicitFactor();
    ExplicitFactor(const ExplicitFactor<T, TABLE>& other);
    ExplicitFactor(const DiscreteSpace&, const T& = static_cast<T>(1));
    template<class IndexIterator>
        ExplicitFactor(const DiscreteSpace&, IndexIterator, IndexIterator, const T& = static_cast<T>(1));
    template<class IndexIterator, class ValueIterator>
        ExplicitFactor(const DiscreteSpace&, IndexIterator, IndexIterator, ValueIterator, ValueIterator);
    ~ExplicitFactor();

    void assign(const DiscreteSpace&, const T& = static_cast<T>(1));
    template<class IndexIterator>
        void assign(const DiscreteSpace&, IndexIterator, IndexIterator, const T& = static_cast<T>(1));
    template<class IndexIterator, class ValueIterator>
        void assign(const DiscreteSpace&, IndexIterator, IndexIterator, ValueIterator, ValueIterator);

    const space_type& space() const;
    bool isConstant() const;
    size_t numberOfVariables() const;
    size_t variableIndex(const size_t&) const;
    bool dependsOnVariable(const size_t&) const;
    template<class Iterator>
    void variableIndices(Iterator) const;
    variable_index_container_type variableIndices() const;
    template<class IndexIterator>
    bool dependsOnVariables(IndexIterator, IndexIterator) const;

    T& operator()();
    template<class U>
        T& operator()(U);
    T& operator()(const size_t&, const size_t&);
    T& operator()(const size_t&, const size_t&, const size_t&);
    T& operator()(const size_t&, const size_t&, const size_t&, const size_t&);
    const T& operator()() const;
    template<class U>
        const T& operator()(U) const;
    const T& operator()(const size_t&, const size_t&) const;
    const T& operator()(const size_t&, const size_t&, const size_t&) const;
    const T& operator()(const size_t&, const size_t&, const size_t&, const size_t&) const;
    const table_type& table() const;

    template<class IndexIterator, class StateIterator>
        void fixVariables(IndexIterator, IndexIterator, StateIterator);

    template<class Accumulator, class IndexIterator>
        void accumulate(IndexIterator, IndexIterator,
            ExplicitFactor<T, TABLE>&) const; // over some variables
    template<class Accumulator>
        void accumulate(T&, std::vector<state_type>&) const; // over all variables

    template<class Accumulator>
        void accumulate(T&) const; // over all variables

    void normalize();
    void subtractOffset();

    ExplicitFactor<T, TABLE>& operator=(const ExplicitFactor<T, TABLE>& other);
    ExplicitFactor<T, TABLE> operator+() const;
    ExplicitFactor<T, TABLE> operator-() const;
    void operator*=(const T&);
    void operator/=(const T&);
    void operator+=(const T&);
    void operator-=(const T&);
    ExplicitFactor<T, TABLE> operator*(const T&) const;
    ExplicitFactor<T, TABLE> operator/(const T&) const;
    ExplicitFactor<T, TABLE> operator+(const T&) const;
    ExplicitFactor<T, TABLE> operator-(const T&) const;
    ExplicitFactor<T, TABLE> operator*(const ExplicitFactor<T, TABLE>&) const;
    ExplicitFactor<T, TABLE> operator/(const ExplicitFactor<T, TABLE>&) const;
    ExplicitFactor<T, TABLE> operator+(const ExplicitFactor<T, TABLE>&) const;
    ExplicitFactor<T, TABLE> operator-(const ExplicitFactor<T, TABLE>&) const;
    void operator*=(const ExplicitFactor<T, TABLE>&);
    void operator/=(const ExplicitFactor<T, TABLE>&);
    void operator+=(const ExplicitFactor<T, TABLE>&);
    void operator-=(const ExplicitFactor<T, TABLE>&);
    void negate();
    void reciprocal();
    template<class U>
    void pow(const U&);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!! the following methods are specific for explicitfactor !!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    void exp();
    void log();
    void abs();

    template<class IndexIterator>
        void integrate(IndexIterator, IndexIterator);
    void integrate(size_t);
    void integrate();

    template<class IndexIterator>
        void maximize(IndexIterator, IndexIterator);
    void maximize(size_t);
    void maximize();

    template<class IndexIterator>
        void minimize(IndexIterator, IndexIterator);
    void minimize(size_t);
    void minimize();

    bool isSubmodular() const;

    // the following methods are public, for convenience

    template<class UnaryFunction>
        void operate(UnaryFunction);
    template<class BinaryFunction>
        void operate(const T&, BinaryFunction);
    template<class BinaryFunction>
        static void operate(const ExplicitFactor<T, TABLE>&, const ExplicitFactor<T, TABLE>&,
            ExplicitFactor<T, TABLE>&, BinaryFunction);
    template<class BinaryFunction>
        static void operate(const ExplicitFactor<T, TABLE>&, ExplicitFactor<T, TABLE>&,
            BinaryFunction);

private:
    const DiscreteSpace* space_;
    variable_index_type* variableIndices_;
    table_type table_;  
};

template<class T, class TABLE>
ExplicitFactor<T, TABLE>::ExplicitFactor
(
    const ExplicitFactor<T, TABLE>& other
)
:   space_(other.space_),
    variableIndices_(new variable_index_type[other.numberOfVariables()]),
    table_(other.table_)
{
    std::copy(other.variableIndices_, other.variableIndices_+other.numberOfVariables(), variableIndices_);
}

template<class T, class TABLE>
ExplicitFactor<T, TABLE>& 
ExplicitFactor<T, TABLE>::operator=
(
    const ExplicitFactor<T, TABLE>& other
) {
    if(&other == this) {
        //handle self-assignment
        return *this;
    }
    space_ = other.space_;
    table_ = other.table_;
    delete[] variableIndices_;
    variableIndices_ = new variable_index_type[other.numberOfVariables()];
    std::copy(other.variableIndices_, other.variableIndices_+other.numberOfVariables(), variableIndices_);

    return *this;
}

template<class T, class TABLE>
void ExplicitFactor<T, TABLE>::assign
(
    const DiscreteSpace& space,
    const T& constant
)
{
    delete[] variableIndices_;
    variableIndices_ = 0;

    space_ = &space;
    table_ = table_type(constant);
    assert(numberOfVariables() == 0);
}

template<class T, class TABLE>
template<class IndexIterator>
inline void
ExplicitFactor<T, TABLE>::assign
(
    const DiscreteSpace& space,
    IndexIterator begin,
    IndexIterator end,
    const T& constant
)
{
    // before we might throw any exceptions, clean up variableIndices_
    delete[] variableIndices_;
    variableIndices_ = 0;
    space_ = 0;

    Assert(NO_ARG_TEST || begin != end);
    Assert(NO_ARG_TEST || space.isValidIndexSequence(begin, end));
    size_t length = std::distance(begin, end);
    //FIXME optimize this without having to do the copy
    variable_index_container_type numbersOfStates(length);
    space.numbersOfStates(begin, end, numbersOfStates.begin());
    table_ = table_type(numbersOfStates.begin(),
        numbersOfStates.end(), constant);
    space_ = &space;

    variableIndices_ = new variable_index_type[length];
    std::copy(begin, end, variableIndices_);
    assert(numberOfVariables() == length);
}

template<class T, class TABLE>
template<class IndexIterator, class ValueIterator>
inline void
ExplicitFactor<T, TABLE>::assign
(
    const DiscreteSpace& space,
    IndexIterator begin,
    IndexIterator end,
    ValueIterator valueBegin,
    ValueIterator valueEnd
)
{
    delete[] variableIndices_;
    variableIndices_ = 0;

    Assert(NO_ARG_TEST || begin != end);
    Assert(NO_ARG_TEST || space.isValidIndexSequence(begin, end));
    size_t length = std::distance(begin, end);
    std::vector<size_t> numbersOfStates(length);
    space.numbersOfStates(begin, end, numbersOfStates.begin());
    table_ = table_type(marray::SkipInitialization,
        numbersOfStates.begin(), numbersOfStates.end());
    for (size_t j=0; valueBegin != valueEnd; ++j, ++valueBegin) {
        table_(j) = *valueBegin;
    }
    space_ = &space;

    variableIndices_ = new variable_index_type[length];
    std::copy(begin, end, variableIndices_);
    assert(numberOfVariables() == length);
}

template<class T, class TABLE>
ExplicitFactor<T, TABLE>::ExplicitFactor()
: space_(0),
  variableIndices_(0),
  table_(table_type())
{
    assert(numberOfVariables() == 0);
}

template<class T, class TABLE>
ExplicitFactor<T, TABLE>::ExplicitFactor
(
    const DiscreteSpace& space,
    const T& constant
)
:   space_(&space),
    variableIndices_(0),
    table_(table_type(constant))
{
    assert(numberOfVariables() == 0);
}

// Construct an explicit potential with the given discrete space
// and set of variable indices
//
template<class T, class TABLE>
template<class IndexIterator>
ExplicitFactor<T, TABLE>::ExplicitFactor
(
    const DiscreteSpace& space,
    IndexIterator begin,
    IndexIterator end,
    const T& constant
)
:   space_(0),
    variableIndices_(0),
    table_(table_type())
{
    assign(space, begin, end, constant);
}

// Construct an explicit potential with the given discrete space,
// set of variable indices, and values
//
template<class T, class TABLE>
template<class IndexIterator, class ValueIterator>
ExplicitFactor<T, TABLE>::ExplicitFactor
(
    const DiscreteSpace& space,
    IndexIterator begin,
    IndexIterator end,
    ValueIterator valueBegin,
    ValueIterator valueEnd
)
:   space_(0),
    variableIndices_(0),
    table_(table_type())
{
    assign(space, begin, end, valueBegin, valueEnd);
}

template<class T, class TABLE>
ExplicitFactor<T, TABLE>::~ExplicitFactor()
{
    delete[] variableIndices_;
    variableIndices_ = 0;
}

template<class T, class TABLE>
inline const typename ExplicitFactor<T, TABLE>::space_type&
ExplicitFactor<T, TABLE>::space() const
{
    Assert(NO_DEBUG || space_ != 0);
    return *space_;
}

template<class T, class TABLE>
inline bool
ExplicitFactor<T, TABLE>::isConstant() const
{
    Assert(NO_DEBUG || space_ != 0);
    return(numberOfVariables() == 0);
}

template<class T, class TABLE>
inline size_t
ExplicitFactor<T, TABLE>::numberOfVariables() const
{
    //FIXME
    //Assert(NO_DEBUG || space_ != 0);
    return table_.size() == 0 ? 0 : table_.dimension();
}

template<class T, class TABLE>
inline size_t
ExplicitFactor<T, TABLE>::variableIndex
(
    const size_t& index
) const
{
    Assert(NO_ARG_TEST || index < numberOfVariables());
    return variableIndices_[index];
}

template<class T, class TABLE>
template<class Iterator>
inline void
ExplicitFactor<T, TABLE>::variableIndices
(
    Iterator out
) const
{
    Assert(NO_DEBUG || space_ != 0);
    for (size_t j=0; j<numberOfVariables(); ++j) {
        *out = variableIndices_[j];
        ++out;
    }
}

template<class T, class TABLE>
inline typename ExplicitFactor<T, TABLE>::variable_index_container_type
ExplicitFactor<T, TABLE>::variableIndices() const
{
    variable_index_container_type v(numberOfVariables());
    std::copy(variableIndices_, variableIndices_+numberOfVariables(), v.begin());
    return v;
}

template<class T, class TABLE>
template<class IndexIterator>
bool
ExplicitFactor<T, TABLE>::dependsOnVariables
(
    IndexIterator begin,
    IndexIterator end
) const
{
    Assert(NO_DEBUG || space_ != 0);
    // copy input sequence to set

    std::vector<variable_index_type> variableIndices(std::distance(begin, end));
    std::copy(begin, end, variableIndices.begin());
    std::sort(variableIndices.begin(), variableIndices.end());

    // check for dependency
    variable_index_type* thisvi = variableIndices_;
    std::vector<variable_index_type>::const_iterator vi = variableIndices.begin();
    while (vi != variableIndices.end()) {
        if (thisvi == variableIndices_+numberOfVariables()) {
            return false;
        }
        if (*vi < *thisvi) {
            return false;
        }
        else if (*vi == *thisvi) {
            ++vi;
            ++thisvi;
        }
        else {
            ++vi;
        }
    }
    return true;
}

template<class T, class TABLE>
inline bool
ExplicitFactor<T, TABLE>::dependsOnVariable
(
    const size_t& variableIndex
) const
{
    Assert(NO_DEBUG || space_ != 0);
    return std::binary_search(variableIndices_, variableIndices_+numberOfVariables(),
        variableIndex);
}

// access to constant
template<class T, class TABLE>
inline T&
ExplicitFactor<T, TABLE>::operator()()
{
    Assert(NO_DEBUG || space_ != 0);
    Assert(NO_ARG_TEST || isConstant());
    return table_(0);
}

// access to constant
template<class T, class TABLE>
inline const T&
ExplicitFactor<T, TABLE>::operator()() const
{
    Assert(NO_DEBUG || space_ != 0);
    Assert(NO_ARG_TEST || isConstant());
    return table_(0);
}

// u can either be an index or an iterator to the beginning of
// a sequence of variable states
template<class T, class TABLE>
template<class U>
inline T&
ExplicitFactor<T, TABLE>::operator()
(
    U u
)
{
    Assert(NO_DEBUG || space_ != 0);
    return table_(u);
}

// u can either be an index or an iterator to the beginning of
// a sequence of variable states
template<class T, class TABLE>
template<class U>
inline const T&
ExplicitFactor<T, TABLE>::operator()
(
    U u
) const
{
    Assert(NO_DEBUG || space_ != 0);
    return table_(u);
}

template<class T, class TABLE>
inline T&
ExplicitFactor<T, TABLE>::operator()
(
    const size_t& index0,
    const size_t& index1
)
{
    return table_(index0, index1);
}

template<class T, class TABLE>
inline T&
ExplicitFactor<T, TABLE>::operator()
(
    const size_t& index0,
    const size_t& index1,
    const size_t& index2
)
{
    return table_(index0, index1, index2);
}

template<class T, class TABLE>
inline T&
ExplicitFactor<T, TABLE>::operator()
(
    const size_t& index0,
    const size_t& index1,
    const size_t& index2,
    const size_t& index3
)
{
    return table_(index0, index1, index2, index3);
}

template<class T, class TABLE>
inline const T&
ExplicitFactor<T, TABLE>::operator()
(
    const size_t& index0,
    const size_t& index1,
    const size_t& index2,
    const size_t& index3
) const
{
    return table_(index0, index1, index2, index3);
}

template<class T, class TABLE>
inline const T&
ExplicitFactor<T, TABLE>::operator()
(
    const size_t& index0,
    const size_t& index1,
    const size_t& index2
) const
{
    return table_(index0, index1, index2);
}

template<class T, class TABLE>
inline const T&
ExplicitFactor<T, TABLE>::operator()
(
    const size_t& index0,
    const size_t& index1
) const
{
    return table_(index0, index1);
}

template<class T, class TABLE>
inline const typename ExplicitFactor<T, TABLE>::table_type&
ExplicitFactor<T, TABLE>::table() const
{
    return table_;
}

template<class T, class TABLE>
template<class IndexIterator, class StateIterator>
void
ExplicitFactor<T, TABLE>::fixVariables
(
    IndexIterator begin,
    IndexIterator end,
    StateIterator state
)
{
    Assert(NO_DEBUG || space_ != 0);
    if (!isConstant()) {
        std::vector<bool> fixIndices(numberOfVariables(), false);
        size_t numberOfFixedVariables = 0;
        std::vector<size_t> offset(numberOfVariables());
        std::vector<size_t> shape(table_.dimension());
        for (size_t j=0; j<table_.dimension(); ++j) {
            shape[j] = table_.shape(j);
        }
        size_t vi=0;
        while (vi < numberOfVariables() && begin != end) {
            if (variableIndices_[vi] < *begin) {
                ++vi;
            }
            else if (variableIndices_[vi] > *begin) {
                ++begin;
                ++state;
            }
            else { // ==
                Assert(NO_ARG_TEST || *state < space_->numberOfStates(variableIndices_[vi]) );
                offset[vi] = *state;
                shape[vi] = 1;
                fixIndices[vi] = true;
                ++numberOfFixedVariables;
                ++vi;
                ++begin;
                ++state;
            }
        }
        size_t numberOfNewVariables =
            numberOfVariables() - numberOfFixedVariables;
        if (numberOfNewVariables == 0) {
            assign(*space_, table_(offset.begin()));
        }
        else {
            // make newVariableIndices explicit:
            variable_index_container_type newVariableIndices(numberOfNewVariables);
            variable_index_container_type oldDimensionOfNewDimension(numberOfNewVariables);
            {
                size_t k = 0;
                for(size_t j=0; j<numberOfVariables(); ++j) {
                    if(!fixIndices[j]) {
                        oldDimensionOfNewDimension[k] = j;
                        newVariableIndices[k] = variableIndices_[j];
                        ++k;
                    }
                }
            }
            // make newNumbersOfStates explicit:
            std::vector<size_t> newNumbersOfStates(numberOfNewVariables);
            for (size_t j=0; j<numberOfNewVariables; ++j) {
                newNumbersOfStates[j] = space_->numberOfStates(newVariableIndices[j]);
            }
            // create new table:
            table_type newTable(newNumbersOfStates.begin(), newNumbersOfStates.end(), size_t(1));
            // copy
            std::vector<size_t> coordinateOld = offset;
            std::vector<size_t> coordinateNew(newNumbersOfStates.size());
            for(;;) {
                // determine old coordinate based on new coordinate
                for(size_t d=0; d<newNumbersOfStates.size(); ++d) {
                    coordinateOld[oldDimensionOfNewDimension[d]] = coordinateNew[d];
                }
                // copy from old to new
                newTable(coordinateNew.begin()) = table_(coordinateOld.begin());
                // increment new coordinate
                size_t d = 0;
                while(d < newNumbersOfStates.size()) {
                    if(coordinateNew[d]+1 < newNumbersOfStates[d]) {
                        ++coordinateNew[d];
                        break;
                    }
                    else {
                        coordinateNew[d] = 0;
                        ++d;
                    }
                }
                if(d == newNumbersOfStates.size()) {
                    break;
                }
            }
            // assign to attributes
            table_ = newTable; //first copy the table, this also sets the correct
            //number of variables
            //now set the new variables indices
            delete[] variableIndices_;
            variableIndices_ = 0;

            variableIndices_ = new variable_index_type[newVariableIndices.size()];
            std::copy(newVariableIndices.begin(), newVariableIndices.end(), variableIndices_);
            assert(numberOfVariables() == newVariableIndices.size());
        }
    }
}

/// executes x = f(x) on every entry x
template<class T, class TABLE>
template<class UnaryFunction>
inline void
ExplicitFactor<T, TABLE>::operate
(
    UnaryFunction f
)
{
    Assert(NO_DEBUG || space_ != 0);
    for (typename table_type::iterator it = table_.begin(); it!=table_.end(); ++it) {
        *it = f(*it);
    }
}

/// executes x = f(x,v) on every entry x
template<class T, class TABLE>
template<class BinaryFunction>
inline void
ExplicitFactor<T, TABLE>::operate
(
    const T& v,
    BinaryFunction f
)
{
    Assert(NO_DEBUG || space_ != 0);
    for (typename table_type::iterator it = table_.begin(); it!=table_.end(); ++it) {
        *it = f(*it, v);
    }
}

/// executes x3 = f(x1,x2) on any relevant pair (x1,x2) of entries
/// of the tables of p1 and p2
template<class T, class TABLE>
template<class BinaryFunction>
void
ExplicitFactor<T, TABLE>::operate //inplace out (op)= in
(
    const ExplicitFactor<T, TABLE>& in,
    ExplicitFactor<T, TABLE>& out,
    BinaryFunction f
)
{
    Assert(NO_ARG_TEST || (in.space_ != 0 && out.space_ != 0 && &(in.space()) == &(out.space()) ) );
    if (in.isConstant() && out.isConstant()) {
        out.operate(in(), f);
    }
    else {
        {
            std::vector<variable_index_type> variableIndicesOut(in.numberOfVariables()+out.numberOfVariables());
            std::copy(in.variableIndices_,
                in.variableIndices_+in.numberOfVariables(),
                variableIndicesOut.begin());
            std::copy(out.variableIndices_,
                out.variableIndices_+out.numberOfVariables(),
                variableIndicesOut.begin()+in.numberOfVariables());
            std::sort(variableIndicesOut.begin(), variableIndicesOut.end());
            std::vector<variable_index_type>::iterator uniqueEnd =
                std::unique(variableIndicesOut.begin(), variableIndicesOut.end());
            const size_t variableIndicesOutSize = std::distance(variableIndicesOut.begin(),uniqueEnd);

            if (out.space_ == 0 || out.numberOfVariables()
                != variableIndicesOutSize
                ) {
                    out.assign(in.space(), variableIndicesOut.begin(), uniqueEnd);
            }
            //variableIndices runs out of scope
        }

        if (in.numberOfVariables() == out.numberOfVariables()) {
            typename table_type::iterator outit = out.table_.begin();
            typename table_type::const_iterator init = in.table_.begin();
            while (outit != out.table_.end()) {
                *outit = f(*outit, *init); // operate
                ++outit;
                ++init;
            }
        }
        else if (out.numberOfVariables()==2 && in.numberOfVariables()==1) {
            if (in.variableIndex(0) == out.variableIndex(0)) {
                for (size_t i=0; i<out.space().numberOfStates(out.variableIndex(0)); ++i) {
                    for (size_t j=0; j<out.space().numberOfStates(out.variableIndex(1)); ++j) {
                        out(i,j) = f(out(i,j), in(i));
                    }
                }
            } else {
                for (size_t i=0; i<out.space().numberOfStates(out.variableIndex(0)); ++i) {
                    for (size_t j=0; j<out.space().numberOfStates(out.variableIndex(1)); ++j) {
                        out(i,j) = f(out(i,j), in(j));
                    }
                }
            }
        }
        else {
            variable_index_container_type correspIndicesIn(in.numberOfVariables());
            variable_index_container_type newOutVariableIndices(out.numberOfVariables());
            std::copy(out.variableIndices_, out.variableIndices_+out.numberOfVariables(), newOutVariableIndices.begin());
            for (size_t j=0; j<in.numberOfVariables(); ++j) {
                for (size_t i=0; i<out.numberOfVariables(); ++i) {
                    if (in.variableIndices_[j] == out.variableIndices_[i]) {
                        correspIndicesIn[j] = i;
                        break;
                    }
                }
                assert(out.variableIndices_[correspIndicesIn[j]] == in.variableIndices_[j]);
            }
            std::vector<size_t> currentStateIn(in.numberOfVariables());
            std::vector<size_t> currentStateOut(out.table_.dimension());
            typename table_type::iterator outit = out.table_.begin();
            while (outit != out.table_.end()) {
                out.table_.indexToCoordinates(std::distance(out.table_.begin(), outit), currentStateOut.begin());
                for (size_t j=0; j<in.numberOfVariables(); ++j) {
                    currentStateIn[j] = currentStateOut[correspIndicesIn[j]];
                }
                *outit = f(*outit, in(currentStateIn.begin())); // operate
                ++outit;
            }
        }
    }
}

/// executes x2 = f(x1,x2) on any relevant pair (x1,x2) of entries of the tables
/// of p1 and p2. x2 have to cover x1 because operation is INPLACE
template<class T, class TABLE>
template<class BinaryFunction>
void
ExplicitFactor<T, TABLE>::operate
(
    const ExplicitFactor<T, TABLE>& in1,
    const ExplicitFactor<T, TABLE>& in2,
    ExplicitFactor<T, TABLE>& out,
    BinaryFunction f
)
{
    Assert(NO_ARG_TEST || (in1.space_ != 0 && in2.space_ != 0 && &(in1.space()) == &(in2.space()) ) );
    if (in1.isConstant()) {
        out = in2;
        out.operate(in1(), f);
    }
    else {
        if (in2.isConstant()) {
            out = in1;
            out.operate(in2(), f);
        }
        else {//in1  and in2 not const
            std::vector<variable_index_type> variableIndicesOut(in1.numberOfVariables()+in2.numberOfVariables());
            std::copy(in1.variableIndices_,
                in1.variableIndices_+in1.numberOfVariables(),
                variableIndicesOut.begin());
            std::copy(in2.variableIndices_,
                in2.variableIndices_+in2.numberOfVariables(),
                variableIndicesOut.begin()+in1.numberOfVariables());
            std::sort(variableIndicesOut.begin(), variableIndicesOut.end());
            std::vector<variable_index_type>::iterator uniqueEnd =
                std::unique(variableIndicesOut.begin(), variableIndicesOut.end());
            const size_t variableIndicesOutSize = std::distance(variableIndicesOut.begin(), uniqueEnd);

            //Resize OutFactor if neccesary
            if (out.space_ == 0  || out.numberOfVariables()
                != variableIndicesOutSize) {
                    out.assign(in1.space(), variableIndicesOut.begin(), uniqueEnd);
            }

            //If factors use equal variables use efficient loop
            if (variableIndicesOutSize == in1.numberOfVariables()
                && in1.numberOfVariables() == in2.numberOfVariables()
                ) {
                    typename table_type::iterator outit = out.table_.begin();
                    typename table_type::const_iterator in1it = in1.table_.begin();
                    typename table_type::const_iterator in2it = in2.table_.begin();
                    while (outit != out.table_.end()) {
                        *outit = f(*in1it,*in2it); // operate
                        ++outit;
                        ++in1it;
                        ++in2it;
                    }
                    return;
            }

            variable_index_container_type correspIndices1(in1.numberOfVariables());
            // variableIndices1[j] equals vairableIndicesP[correspIndices1[j]]
            variable_index_container_type correspIndices2(in2.numberOfVariables());
            // variableIndices2[j] equals vairableIndicesP[correspIndices1[j]]
            size_t j1=0; // indices for array p1.variableIndices
            size_t j2=0; // indices for array p2.variableIndices
            size_t count=0; // counter of indices of the product
            // determine variableIndicesOut, correspIndices1, correspIndices2
            for (;;) {
                if (j1 == in1.numberOfVariables()) {
                    while (j2 != in2.numberOfVariables()) {
                        //variableIndicesOut.insert(in2.variableIndices_[j2]);
                        correspIndices2[j2] = count;
                        ++j2;
                        ++count;
                    }
                    break;
                }
                else {
                    if (j2 == in2.numberOfVariables()) {
                        while (j1 != in1.numberOfVariables()) {
                            //variableIndicesOut.insert(in1.variableIndices_[j1]);
                            correspIndices1[j1] = count;
                            ++j1;
                            ++count;
                        }
                        break;
                    }
                    else {
                        if (in1.variableIndices_[j1] < in2.variableIndices_[j2]) {
                            //variableIndicesOut.insert(in1.variableIndices_[j1]);
                            correspIndices1[j1] = count;
                            ++j1;
                            ++count;
                        }
                        else if (in1.variableIndices_[j1] > in2.variableIndices_[j2]) {
                            //variableIndicesOut.insert(in2.variableIndices_[j2]);
                            correspIndices2[j2] = count;
                            ++j2;
                            ++count;
                        }
                        else { // ==
                            //variableIndicesOut.insert(in1.variableIndices_[j1]);
                            correspIndices1[j1] = count;
                            correspIndices2[j2] = count;
                            ++j1;
                            ++j2;
                            ++count;
                        }
                    }
                }
            }
            std::vector<size_t> currentStateIn1(in1.numberOfVariables());
            std::vector<size_t> currentStateIn2(in2.numberOfVariables());
            typename table_type::iterator outit = out.table_.begin();
            while (outit != out.table_.end()) {
                std::vector<size_t> currentStateOut(out.table_.dimension());
                out.table_.indexToCoordinates(std::distance(out.table_.begin(), outit), currentStateOut.begin());
                for (size_t j=0; j<in1.numberOfVariables(); ++j) {
                    currentStateIn1[j] = currentStateOut[correspIndices1[j]];
                }
                for (size_t j=0; j<in2.numberOfVariables(); ++j) {
                    currentStateIn2[j] = currentStateOut[correspIndices2[j]];
                }
                *outit = f(in1(currentStateIn1.begin()),
                    in2(currentStateIn2.begin())); // operate
                ++outit;
            }
        }
    }
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::negate()
{
    operate(std::negate<T>());
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::reciprocal()
{
    operate(Reciprocal<T>());
}

template<class T, class TABLE>
template<class Accumulator>
void ExplicitFactor<T, TABLE>::accumulate
(
    T& value
) const
{
    typedef typename table_type::const_iterator Ti;
    Assert(NO_DEBUG || space_ != 0);
    Accumulation<value_type,state_type, Accumulator> acc;
    if (isConstant()) {
        acc(table_(0));
    }
    else {
        for (Ti ti = table_.begin(); ti!=table_.end(); ++ti) {
            acc(*ti);
        }
    }
    value = acc.value();
}

template<class T, class TABLE>
template<class Accumulator>
void ExplicitFactor<T, TABLE>::accumulate
(
    T& value,
    std::vector<state_type>& state
) const
{
    Assert(NO_DEBUG || space_ != 0);
    Accumulation<value_type, state_type, Accumulator> acc;
    if (isConstant()) {
        std::vector<state_type> emptyState;
        acc(table_(0), emptyState);
    }
    else {
        for (typename table_type::const_iterator ti = table_.begin(); ti!=table_.end(); ++ti) {
            std::vector<size_t> seq(table_.dimension());
            const size_t d = std::distance(table_.begin(), ti);
            table_.indexToCoordinates(d, seq.begin());
            acc(*ti, seq);
        }
    }
    value = acc.value();
    acc.state(state);
}

template<class T, class TABLE>
template<class Accumulator, class IndexIterator>
void ExplicitFactor<T, TABLE>::accumulate
(
    IndexIterator begin,
    IndexIterator end,
    ExplicitFactor<T, TABLE>& out
) const
{
    Assert(NO_DEBUG || space_ != 0);

    // copy passed indices,
    // sort them, and do a set intersection with variableIndices_.
    // Afterwards we need to remove possible duplicate elements
    std::vector<variable_index_type> accumulationIndicesSet(std::distance(begin, end));
    std::copy(begin, end, accumulationIndicesSet.begin());
    std::sort(accumulationIndicesSet.begin(), accumulationIndicesSet.end());
    variable_index_container_type accumulationIndices(accumulationIndicesSet.size());
    variable_index_container_type::iterator ait = accumulationIndices.begin();
    ait = std::set_intersection(accumulationIndicesSet.begin(), accumulationIndicesSet.end(),
        variableIndices_, variableIndices_+numberOfVariables(), ait);
    ait = std::unique(accumulationIndices.begin(), ait);
    size_t accumulationIndicesSize = std::distance(accumulationIndices.begin(), ait);

    // for each dimension of the table, identify whether or not to
    // accumulate over this dimension
    std::vector<bool> accumulateDimension(numberOfVariables(), false);
    {
        size_t k=0;
        for(size_t j=0; j<accumulationIndicesSize; ++j) {
            while(variableIndices_[k] != accumulationIndices[j]) {
                ++k;
            }
            accumulateDimension[k] = true;
        }
    }

    std::vector<size_t> newNumbersOfStates(numberOfVariables() - accumulationIndicesSize);
    variable_index_container_type accVariableIndices(accumulationIndicesSize);
    // indices of variables accumulated over
    variable_index_container_type notAccVariableIndices(newNumbersOfStates.size());
    // indices of variables not accumulated over
    size_t k = 0;
    size_t m = 0;
    for (size_t j=0; j<numberOfVariables(); ++j) {
        if (accumulateDimension[j]) {
            accVariableIndices[m] = j;
            ++m;
        }
        else {
            notAccVariableIndices[k] = variableIndices_[j];
            newNumbersOfStates[k] = space_->numberOfStates(variableIndices_[j]);
            ++k;
        }
    }

    // initialize output, if necessary
    {
        bool initOutput = false;
        if (out.space_ != space_ || k != out.numberOfVariables()) {
            initOutput = true;
        }
        else {
            for (size_t j=0; j<k; ++j) {
                if (notAccVariableIndices[j] != out.variableIndex(j)) {
                    initOutput = true;
                    break;
                }
            }
        }
        if (initOutput) {
            if (notAccVariableIndices.size() == 0) {
                out = ExplicitFactor<T, TABLE>(space());
            }
            else {
                out = ExplicitFactor<T, TABLE>(space(), notAccVariableIndices.begin(), notAccVariableIndices.end());
            }
        }
    }

    // accumulation
    if (notAccVariableIndices.size() == 0) {
        std::vector<size_t> states;
        accumulate<Accumulator>(out(0), states);
    } else if (notAccVariableIndices.size()==1 && accVariableIndices.size()==1 ) {
        if (variableIndex(0) == out.variableIndex(0)) {
            for (size_t i=0; i<space_->numberOfStates(variableIndex(0)); ++i) {
                Accumulator::neutral(out(i));
                for (size_t j=0; j<space_->numberOfStates(variableIndex(1)); ++j) {
                    Accumulator::op(out(i),table_(i,j),out(i));
                }
            }
        } else {
            for (size_t i=0; i<space_->numberOfStates(variableIndex(1)); ++i) {
                Accumulator::neutral(out(i));
                for (size_t j=0; j<space_->numberOfStates(variableIndex(0)); ++j) {
                    Accumulator::op(out(i),table_(j,i),out(i));
                }
            }
        }
    } else {
        marray::Marray<Accumulation<value_type, state_type, Accumulator> >
            acc(newNumbersOfStates.begin(), newNumbersOfStates.end());
        typename marray::Marray<Accumulation<value_type, state_type, Accumulator> >::iterator
            ai = acc.begin();
        std::vector<size_t> offset(numberOfVariables());
        std::vector<size_t> shape(numberOfVariables());
        // iterate over entries of out
        size_t i = 0;
        while (ai.hasMore()) {
            // determine start and end coordinates in table:
            size_t k = 0;
            for (size_t j=0; j<numberOfVariables(); ++j) {
                if (accumulateDimension[j]) {
                    offset[j] = 0;
                    shape[j] = space_->numberOfStates(variableIndices_[j]);
                }
                else {
                    std::vector<size_t> coordinate(acc.dimension());
                    ai.coordinate(coordinate.begin());
                    offset[j] = coordinate[k];
                    shape[j] = 1;
                    ++k;
                }
            }

            // copy
            std::vector<size_t> coordinateTable = offset;
            std::vector<size_t> coordinateAcc(accVariableIndices.size());
            for(;;) {
                // determine coordinate in table_ based on coordinateAcc
                for(size_t d=0; d<coordinateAcc.size(); ++d) {
                    coordinateTable[accVariableIndices[d]] 
                    = coordinateAcc[d]+offset[accVariableIndices[d]] ;
                }
                (*ai)(table_(coordinateTable.begin()), coordinateAcc); 
                // increment coordinateAcc
                size_t d = 0;
                while(d < coordinateAcc.size()) {
                    if(coordinateAcc[d]+1 < shape[accVariableIndices[d]]) {
                        ++coordinateAcc[d];
                        break;
                    }
                    else {
                        coordinateAcc[d] = 0;
                        ++d;
                    }
                }
                if(d == coordinateAcc.size()) {
                    break;
                }
            }
            out(i) = (*ai).value();
            ++ai;
            ++i;
        }
    }
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::normalize()
{
    Assert(NO_DEBUG || space_ != 0);
    if (isConstant()) {
        table_(0) = static_cast<T>(1);
    }
    else {
        T v;
        std::vector<state_type> states;
        accumulate<Integrator>(v, states);
        *this /= v;
    }
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::subtractOffset()
{
    Assert(NO_DEBUG || space_ != 0);
    if (isConstant()) {
        table_(0) = static_cast<T>(0);
    }
    else {
        T v;
        std::vector<state_type>  states;
        accumulate<Minimizer>(v,states);
        *this -= v;
    }
}

template<class T, class TABLE>
template<class IndexIterator>
inline void ExplicitFactor<T, TABLE>::maximize
(
    IndexIterator begin,
    IndexIterator end
)
{
    ExplicitFactor<T, TABLE> temp;
    accumulate<Maximizer, IndexIterator>(begin, end, temp);
    (*this) = temp;
}

template<class T, class TABLE>
inline void ExplicitFactor<T, TABLE>::maximize
(
    size_t variableIndex
)
{
    maximize(&variableIndex, (&variableIndex)+1);
}

template<class T, class TABLE>
inline void ExplicitFactor<T, TABLE>::maximize()
{
    maximize(variableIndices_, variableIndices_+numberOfVariables());
}

template<class T, class TABLE>
template<class IndexIterator>
inline void ExplicitFactor<T, TABLE>::minimize
(
    IndexIterator begin,
    IndexIterator end
)
{
    ExplicitFactor<T, TABLE> temp;
    accumulate<Minimizer, IndexIterator>(begin, end, temp);
    (*this) = temp;
}

template<class T, class TABLE>
inline void ExplicitFactor<T, TABLE>::minimize
(
    size_t variableIndex
)
{
    minimize(&variableIndex, (&variableIndex)+1);
}

template<class T, class TABLE>
inline void ExplicitFactor<T, TABLE>::minimize()
{
    minimize(variableIndices_, variableIndices_+numberOfVariables());
}

template<class T, class TABLE>
template<class IndexIterator>
inline void ExplicitFactor<T, TABLE>::integrate
(
    IndexIterator begin,
    IndexIterator end
)
{
    ExplicitFactor<T, TABLE> temp;
    accumulate<Integrator, IndexIterator>(begin, end, temp);
    (*this) = temp;
}

template<class T, class TABLE>
inline void ExplicitFactor<T, TABLE>::integrate
(
    size_t variableIndex
)
{
    maximize(&variableIndex, (&variableIndex)+1);
}

template<class T, class TABLE>
inline void ExplicitFactor<T, TABLE>::integrate()
{
    maximize(variableIndices_, variableIndices_+numberOfVariables());
}

template<class T, class TABLE>
bool
ExplicitFactor<T, TABLE>::isSubmodular() const
{
    if (!NO_DEBUG) {
        if (numberOfVariables() > 2) {
            throw std::runtime_error("function isSubmodular() is only implemented for factors of order 2.");
        }
        for (size_t j=0; j<numberOfVariables(); ++j) {
            if (space_->numberOfStates(variableIndices_[j]) != 2) {
                throw std::runtime_error("function isSubmodular() is only implemented for binary variables.");
            }
        }
    }
    if (numberOfVariables() < 2) {
        return true;
    }
    else {
        return (*this)(0,0) + (*this)(1,1) <= (*this)(0,1) + (*this)(1,0);
    }
}

// arithmetic operators

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
operator+
(
    const T& v,
    const ExplicitFactor<T, TABLE>& f
)
{
    return f + v;
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
operator-
(
    const T& v,
    const ExplicitFactor<T, TABLE>& f
)
{
    return v + (-f);
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
operator*
(
    const T& v,
    const ExplicitFactor<T, TABLE>& f
)
{
    return f * v;
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
operator/
(
    const T& v,
    const ExplicitFactor<T, TABLE>& f
)
{
    ExplicitFactor<T, TABLE> tmp = f;
    tmp.reciprocal();
    return v * tmp;
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
ExplicitFactor<T, TABLE>::operator+() const
{
    return *this;
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
ExplicitFactor<T, TABLE>::operator-() const
{
    ExplicitFactor<T, TABLE> out = *this;
    out.operate(std::negate<T>());
    return out;
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::operator*=(const T& v)
{
    operate(v, std::multiplies<T>());
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::operator/=(const T& v)
{
    operate(v, std::divides<T>());
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::operator+=(const T& v)
{
    operate(v, std::plus<T>());
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::operator-=(const T& v)
{
    operate(v, std::minus<T>());
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::operator*=(const ExplicitFactor<T, TABLE>& p)
{
    operate(p, *this, std::multiplies<T>());
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::operator/=(const ExplicitFactor<T, TABLE>& p)
{
    operate(p, *this, std::divides<T>());
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::operator+=(const ExplicitFactor<T, TABLE>& p)
{
    operate(p, *this, std::plus<T>());
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::operator-=(const ExplicitFactor<T, TABLE>& p)
{
    operate(p, *this, std::minus<T>());
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
ExplicitFactor<T, TABLE>::operator*(const T& v) const
{
    ExplicitFactor<T, TABLE> out = *this;
    out *= v;
    return out;
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
ExplicitFactor<T, TABLE>::operator/(const T& v) const
{
    ExplicitFactor<T, TABLE> out = *this;
    out /= v;
    return out;
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
ExplicitFactor<T, TABLE>::operator+(const T& v) const
{
    ExplicitFactor<T, TABLE> out = *this;
    out += v;
    return out;
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
ExplicitFactor<T, TABLE>::operator-(const T& v) const
{
    ExplicitFactor<T, TABLE> out = *this;
    out -= v;
    return out;
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
ExplicitFactor<T, TABLE>::operator*(const ExplicitFactor<T, TABLE>& p) const
{
    ExplicitFactor<T, TABLE> out;
    operate(*this, p, out, std::multiplies<T>());
    return out;
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
ExplicitFactor<T, TABLE>::operator/(const ExplicitFactor<T, TABLE>& p) const
{
    ExplicitFactor<T, TABLE> out;
    operate(*this, p, out, std::divides<T>());
    return out;
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
ExplicitFactor<T, TABLE>::operator+(const ExplicitFactor<T, TABLE>& p) const
{
    ExplicitFactor<T, TABLE> out;
    operate(*this, p, out, std::plus<T>());
    return out;
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
ExplicitFactor<T, TABLE>::operator-
(
    const ExplicitFactor<T, TABLE>& p
) const
{
    ExplicitFactor<T, TABLE> out;
    operate(*this, p, out, std::minus<T>());
    return out;
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::abs()
{
    for (size_t j=0; j<table_.size(); ++j) {
        table_(j) = std::abs(table_(j));
    }
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
abs
(
    const ExplicitFactor<T, TABLE>& f
)
{
    ExplicitFactor<T, TABLE> tmp = f;
    tmp.abs();
    return tmp;
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::exp()
{
    for (size_t j=0; j<table_.size(); ++j) {
        table_(j) = std::exp(table_(j));
    }
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
exp
(
    const ExplicitFactor<T, TABLE>& f
)
{
    ExplicitFactor<T, TABLE> tmp = f;
    tmp.exp();
    return tmp;
}

template<class T, class TABLE>
inline void
ExplicitFactor<T, TABLE>::log()
{
    for (size_t j=0; j<table_.size(); ++j) {
        table_(j) = std::log(table_(j));
    }
}

template<class T, class TABLE>
inline ExplicitFactor<T, TABLE>
log
(
    const ExplicitFactor<T, TABLE>& f
)
{
    ExplicitFactor<T, TABLE> tmp = f;
    tmp.log();
    return tmp;
}

template<class T, class TABLE>
template<class U>
inline void
ExplicitFactor<T, TABLE>::pow
(
    const U& u
)
{
    for (size_t j=0; j<table_.size(); ++j) {
        table_(j) = std::pow(table_(j), u);
    }
}

template<class T, class TABLE, class U>
inline ExplicitFactor<T, TABLE>
pow
(
    const ExplicitFactor<T, TABLE> f,
    const U& u
)
{
    ExplicitFactor<T, TABLE> tmp = f;
    tmp.pow(u);
    return tmp;
}

} // namespace opengm

#endif // #ifndef OPENGM_EXPLICITFACTOR_HXX
