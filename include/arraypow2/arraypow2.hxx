/// ArrayPow2. Copyright (c) 2010 by Thorben Kroeger.
///
/// This software was developed by Thorben Kroeger.
/// Enquiries shall be directed to: thorben.kroeger@iwr.uni-heidelberg.de
///
/// All advertising materials mentioning features or use of this software must
/// display the following acknowledgement: ``This product includes ArrayPow2
/// developed by Thorben Kroeger. Please direct enquiries concerning ArrayPow2
/// to thorben.kroeger@iwr.uni-heidelberg.de''.
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
///   display the following acknowledgement: ``This product includes ArrayPow2
///   developed by Thorben Kroeger. Please direct enquiries concerning ArrayPow2
///   to thorben.kroeger@iwr.uni-heidelberg.de''.
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
#ifndef ARRAYPOW2_HXX
#define ARRAYPOW2_HXX

#include <cassert>
#include <limits>
#include <ostream>
#include <marray/marray.hxx>

inline size_t pow2(size_t in)
{ return static_cast<size_t>(1) << in; }

template<bool isIntegral> struct AccessOperatorHelper;
template<class T> class ArrayPow2;
template<class T>
std::ostream& operator<<(std::ostream& o, const ArrayPow2<T>& array);

/// Multi-dimensional arrays whose shape is two for each dimension
///
/// In the case of binary variables, the value table of a function
/// can be represented by a multidimensional table of type T which
/// has two entries in each dimension.
/// ArrayPow2 is heavily optimized for this case, using the least amount
/// of memory possible and utilizing fast binary operations throughout
/// the code.
template<class T>
class ArrayPow2  {
public:
    typedef const T* const_iterator;
    typedef T* iterator;
    typedef T value_type;
        
    friend struct AccessOperatorHelper<true>;
    friend struct AccessOperatorHelper<false>;
    friend std::ostream& operator<<<>(std::ostream& o, const ArrayPow2<T>& array);
    
    inline ArrayPow2() 
    : dimension_(0), data_(0) 
    {}
    
    inline ArrayPow2(const ArrayPow2<T>& other)
    : dimension_(other.dimension_), data_(0)
    {
        if(other.size() != 0) {
            assert(other.data_ != 0);
            assert(other.size() == pow2(other.dimension_));
            data_ = new T[other.size()];
            memcpy(data_, other.data_, sizeof(T)*other.size());
        }
        else {
            assert(other.dimension_ == 0);
        }
    }
    
    inline ~ArrayPow2() { 
        delete[] data_; 
    }
    
    template<class Iter>
    inline ArrayPow2(const marray::InitializationSkipping, Iter begin, Iter end)
    : dimension_(std::distance(begin, end)), data_(0) {
        data_ = new T[ pow2(dimension_) ];
    }
    
    template<class Iter>
    inline ArrayPow2(Iter begin, Iter end, const T& init = T())
    : dimension_(std::distance(begin, end)), data_(0) {
        data_ = new T[ pow2(dimension_) ];
        for(size_t j=0; j<size(); ++j) {
            data_[j] = init;
        }
    }

    // initialize as constant
    inline ArrayPow2(const value_type& value)
    : dimension_(0), data_(0) {
        data_ = new T[1];
        data_[0] = value;
    }

    size_t size() const {
        if(data_ == 0) { 
            // array is un-initialized
            assert(dimension_ == 0);
            return 0;
        }
        else {
            return pow2(dimension_);
            // note: if the array is a constant and thus,
            // dimension_ == 0, the desired size, 1, is returned
        }
    }
    
    inline ArrayPow2<T>& operator=(const ArrayPow2<T>& other) {
        if(&other == this) {
            //handle self-assignment
            return *this;
        }
        else {
            if(other.data_ == 0) {
                delete[] data_;
                data_ = 0;
                assert(other.dimension_ == 0);
                dimension_ = 0;
            }
            else {
                assert(other.size() > 0 && other.size() == pow2(other.dimension_));
                if(size() != other.size()) {
                    delete[] data_;
                    data_ = new T[other.size()];
                }
                memcpy(data_, other.data_, sizeof(T)*other.size());
                dimension_ = other.dimension_;
            }
            return *this;
        }
    }
    
    inline bool operator==(const ArrayPow2<T>& other) const {
        if(dimension_ != other.dimension_) {
            return false;
        }
        else {
            assert(size() == other.size());
            for(size_t i=0; i<size(); ++i) {
                if(data_[i] != other.data_[i]) {
                    return false;
                }
            }
            return true;
        }
    }
    
    template<class Iter>
    inline void coordinateToIndex(Iter b, size_t& out) const {
        /*
        out = 0;
        Iter begin = b;
        if(dimension_ <= 0) { return; }
        for(size_t i=dimension_; i>=1; --i) {
            out |= ( *begin << (i-1));
            ++begin;
        }
        */
        out = 0;
        if(dimension_ == 0) {
            return;
        }
        else {
            for(size_t j=0; j<dimension_; ++j, ++b) {
                out |= static_cast<size_t>(*b) << j;
            }
        }
    }
    
    template<class Iter>
    inline void indexToCoordinates(const size_t index, Iter begin) const {
        /*
        for(size_t i=dimension_; i>=1; --i, ++begin) {
            *begin = (index & (1 << (i-1))) != 0 ? 1 : 0;
        }
        */
        for(size_t j=0; j<dimension_; ++j, ++begin) {
            *begin = (index & pow2(j)) != 0 ? 1 : 0;
        }
    }

    //unary
    template<class U>
    inline const T& operator()(const U x) const {
        return AccessOperatorHelper<std::numeric_limits<U>::is_integer>::execute(*this, x);
    }
    template<class U>
    inline T& operator()(const U x) {
        return AccessOperatorHelper<std::numeric_limits<U>::is_integer>::execute(*this, x);
    }

    //binary
    inline const T& operator()(const size_t x, const size_t y) const {
        return data_[ (x<<0) | (y<<1) ];
    }
    inline T& operator()(const size_t x, const size_t y) {
        return data_[ (x<<0) | (y<<1) ];
    }
    
    //third order
    inline const T& operator()(const size_t x, const size_t y, const size_t z) const {
        return data_[ (x<<0) | (y<<1) | (z<<2) ];
    }
    inline T& operator()(const size_t x, const size_t y, const size_t z) {
        return data_[ (x<<0) | (y<<1) | (z<<2) ];
    }
    
    //fourth order
    inline const T& operator()(const size_t x, const size_t y, const size_t z, const size_t a) const {
        return data_[ (x<<0) | (y<<1) | (z<<2) | (a<<3) ];
    }
    inline T& operator()(const size_t x, const size_t y, const size_t z, const size_t a) {
        return data_[ (x<<0) | (y<<1) | (z<<2) | (a<<3) ];
    }
    
    inline const T* begin() const {
        return data_;
    }
    inline const T* end() const {
        return data_+size();
    }
    inline T* begin() {
        return data_;
    }
    inline T* end() {
        return data_+size();
    }
    
    inline size_t dimension() const { return static_cast<size_t>(dimension_); }
    inline size_t shape(size_t dimension) const { return 2; }

private:
    unsigned char dimension_;
    T* data_;
};

template<class T>
std::ostream& operator<<(std::ostream& o, const ArrayPow2<T>& array) {
    o << "{" << std::endl;
    for(size_t i=0; i<array.size(); ++i) {
        size_t* c = new size_t[array.dimension()];
        array.indexToCoordinates(i, c);
        o << "\t(";
        for(unsigned char j=0; j<array.dimension_; ++j) {
            o << c[j];
            if(j!=array.dimension_-1) {
                o << ", ";
            }
        }
        o << ") = " << array.data_[i] << std::endl;
        delete[] c;
    }
    o << "}" << std::endl;
    
    return o;
}

template<>
struct AccessOperatorHelper<true>
{
    // access by scalar index
    template<class T>
    inline static T&
    execute(const ArrayPow2<T>& array, const size_t index)
    {
        return array.data_[index];
    }
};

template<>
struct AccessOperatorHelper<false>
{
    // access by iterator
    template<class T, class U>
    inline static T&
    execute(const ArrayPow2<T>& array, const U iterator)
    {
        size_t index = 0;
        array.coordinateToIndex(iterator, index);
        return array.data_[index];
    }
};

#endif //ARRAYPOW2_HXX
