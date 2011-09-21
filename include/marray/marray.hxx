/// \mainpage
/// Marray: Runtime-Flexible Multi-dimensional Views and Arrays in C++.
///
/// \section section_licence Licence
///
/// Copyright (c) 2010 Bjoern Andres.
/// 
/// This software was developed by Bjoern Andres.
/// Enquiries shall be directed to bjoern@andres.sc.
///
/// All advertising materials mentioning features or use of this software must
/// display the following acknowledgement: ``This product includes the Marray 
/// package developed by Bjoern Andres. Please direct enquiries concerning the 
/// Marray package to bjoern@andres.sc''.
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
///   display the following acknowledgement: ``This product includes the Marray 
///   package developed by Bjoern Andres. Please direct enquiries concerning the 
///   Marray package to bjoern@andres.sc''.
/// - The name of the author must not be used to endorse or promote products 
///   derived from this software without specific prior written permission.
///
/// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
/// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
/// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
/// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
/// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
/// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
/// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
/// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
/// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/// 
/// \section section_features Features
/// - Multi-dimensional arrays with runtime-flexible
///   - dimension
///   - shape and size
///   - indexing order (first coordinate major order and last coordinate major order)
///   .
/// - Derived classes for matrices and vectors
/// - Access to entries by
///   - coordinates
///   - scalar indices
///   - STL compliant random access iterators
///   .
/// - Multi-dimensional views on
///   - subsets of multi-dimensional arrays
///   - intervals of memory
///   .
/// 
/// \section section_assertion_testing Assertions and Invariant Testing
/// - Set NO_DEBUG = true to disable general consistency tests.
///   This is the default if NDEBUG is set.
/// - Set NO_ARG_TEST = true to disable function argument tests.
///   This is the default if NDEBUG is set.
/// - Tests can significantly affect runtime performance.
/// 
/// \section section_tested_compilers Tested Compilers
/// - Microsoft Visual Studio 2005 Professional,
///   Windows XP Professional, Service Pack 3, 32 bit
/// - GNU gcc 4.3.2, Linux Debian 4.3.2-1.1, 64 bit
/// 
/// \section section_terminology Terminology
/// - Offset: offset added to a memory address
/// - Index: scalar scan-order index in a (possibly strided) View
/// - Coordinate: vectorial coordinate in a View
/// - Simple array: an array that is unstrided (i.e. the strides
///   equal the shape strides) and has a zero offset (cf. the
///   function testInvariant of View for the formal definition).
///
/// \section section_cpp0x C++0x Extensions
/// - C++0x Extensions are enabled by defining
///   - HAVE_CPP0X_TEMPLATE_TYPEDEFS
///   - HAVE_CPP0X_VARIADIC_TEMPLATES
///   - HAVE_CPP0X_INITIALIZER_LISTS
///   .
/// 

#pragma once
#ifndef MARRAY_HXX
#define MARRAY_HXX

#define MARRAY_COMPATIBILITY 1 // compatiblity with Vigra MultiArray

#include <cassert>
#include <stdexcept> // runtime_error
#include <limits>
#include <string>
#include <sstream>
#include <cstring> // memcpy
#include <iterator> // reverse_iterator, distance
#include <vector>
#include <set>
#include <iostream> // cout

/// Marray namespace.
namespace marray {

// assertion testing

#ifdef NDEBUG
	const bool NO_DEBUG = true; ///< General assertion testing disabled.
	const bool NO_ARG_TEST = true; ///< Argument testing disabled.
#else
	const bool NO_DEBUG = false; ///< General assertion testing enabled.
	const bool NO_ARG_TEST = false; ///< Argument testing enabled.
#endif

/// Assertion testing.
template<class A> inline void Assert(A assertion) {
	if(!assertion) throw std::runtime_error("Assertion failed.");
}

// meta-programming

// \cond suppress doxygen
template <bool PREDICATE, class TRUECASE, class FALSECASE>
	struct IfBool;
template <class TRUECASE, class FALSECASE>
	struct IfBool<true, TRUECASE, FALSECASE>
	{ typedef TRUECASE type; };
template <class TRUECASE, class FALSECASE>
	struct IfBool<false, TRUECASE, FALSECASE>
	{ typedef FALSECASE type; };
// \endcond end suppress doxygen

// namespace types

enum StringStyle {TableStyle, MatrixStyle}; ///< Flag to be used with the member function asString() of View.
enum CoordinateOrder {FirstMajorOrder, LastMajorOrder}; ///< Flag setting the order of coordinate tuples.
struct InitializationSkipping { }; ///< Flag to indicate initialization skipping.

// namespace variables

static const bool Const = true; ///< Flag to be used with the template parameter isConst of View and Iterator.
static const bool Mutable = false; ///< Flag to be used with the template parameter isConst of View and Iterator.
static const CoordinateOrder defaultOrder = LastMajorOrder; ///< Default order of coordinate tuples.
static const InitializationSkipping SkipInitialization = InitializationSkipping(); ///< Flag to indicate initialization skipping.

// namespace functions

template<class ShapeIterator, class StridesIterator>
	void shapeStrides(ShapeIterator, ShapeIterator,
		StridesIterator, const CoordinateOrder& = defaultOrder);

// class template declarations

template<class T> class Vector;
template<class T> class Matrix;
template<class T> class Marray;
template<class T, bool isConst = false> class View;
// \cond suppress doxygen
template<class T, bool isConst> struct AssignmentOperatorHelper;
template<bool isIntegral> struct AccessOperatorHelper;
// \endcond end suppress doxygen
template<class T, bool isConst> class Iterator;

// template typedefs (C++0x)

#ifdef HAVE_CPP0X_TEMPLATE_TYPEDEFS
	template<class T> using ConstView = View<T, true>;
#endif

// classes

/// Interface to an interval of memory that behaves like Marray.
///
/// View can be used to interface data in an interval of memory
/// as if this data was stored in a Marray. Data items can be
/// referenced by coordinates, by scalar indices, as well as by
/// STL compliant random access iterators.
///
/// Notes on arithmetic operators on View:
///
/// For views, the postfix operators ++ and -- cannot be
/// well-defined because the return value of these operators
/// would have to be the view as it is prior to the operator call.
/// However, the data under the view cannot be preserved when 
/// incrementing or decrementing. Thus, the postfix operators ++
/// and -- are not defined for views. Compilers might use the prefix
/// variant implicitly and issue a warning.
/// 
/// The binary operators += -= *= /= are defined for Views
/// which possibly have a different value_type (template
/// parameter T). In constrast, the operators + - * / only work
/// with views having the same value type. This is because
/// the return type of + - * / should not depend on the order of
/// their arguments.
///
/// Notes on C++0x extensions:
///
/// C++0x has support for varidic templates. We use this to allow an
/// arbitrary number of parameters in operator(). If the number of
/// arguments supplied to operator() does not match the dimension, a
/// runtime error is issued.
///
template<class T, bool isConst> 
class View
{
public:
	typedef T value_type;
	typedef typename IfBool<isConst, const T*, T*>::type pointer;
	typedef const T* const_pointer;
	typedef typename IfBool<isConst, const T&, T&>::type reference;
	typedef const T& const_reference;
	typedef Iterator<T, isConst> iterator;
	typedef Iterator<T, true> const_iterator;
	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
	typedef unsigned short dimension_type;
	typedef size_t coordinate_type;
	typedef std::vector<coordinate_type> coordinate_tuple;
	typedef std::vector<size_t> offset_tuple;
	#ifdef MARRAY_COMPATIBILITY
		typedef std::vector<coordinate_type> difference_type; ///< For compatibility with Vigra MultiArray
		typedef coordinate_type difference_type1; ///< For compatibility with Vigra MultiArray
	#endif

	// construction
  void assign();
	template<class ShapeIterator>
		void assign(ShapeIterator, ShapeIterator, pointer,
			const CoordinateOrder& = defaultOrder,
			const CoordinateOrder& = defaultOrder);
	template<class ShapeIterator, class StrideIterator>
		void assign(ShapeIterator, ShapeIterator, StrideIterator,
			pointer, const size_t&, const CoordinateOrder&);
	View();
	View(pointer); 
	View(const View<T, false>&);
	template<class ShapeIterator>
		View(ShapeIterator, ShapeIterator, pointer,
			const CoordinateOrder& = defaultOrder,
			const CoordinateOrder& = defaultOrder);
	template<class ShapeIterator, class StrideIterator>
		View(ShapeIterator, ShapeIterator, StrideIterator,
			pointer, const size_t&, const CoordinateOrder&);
	#ifdef HAVE_CPP0X_INITIALIZER_LISTS
		void assign(std::initializer_list<size_t>, pointer,
			const CoordinateOrder& = defaultOrder,
			const CoordinateOrder& = defaultOrder);
		void assign(std::initializer_list<size_t>, 
			std::initializer_list<size_t>, pointer,
			const size_t&, const CoordinateOrder&);
		View(std::initializer_list<size_t>, pointer,
			const CoordinateOrder& = defaultOrder,
			const CoordinateOrder& = defaultOrder);
		View(std::initializer_list<size_t>, std::initializer_list<size_t>,
			pointer, const size_t&, const CoordinateOrder&);
	#endif
	
	// assignment 
	// it is necessary to define both operators in order to
	// overwrite the standard operator= for each isConst
	View<T, isConst>& operator=(const View<T, true>&); 
	View<T, isConst>& operator=(const View<T, false>&);
	View<T, isConst>& operator=(const T&);

	// query
	const dimension_type& dimension() const;
	const size_t& size() const;
	const coordinate_type& shape(const dimension_type&) const;
	const coordinate_tuple& shape() const;
	const size_t& strides(const dimension_type&) const;
	const offset_tuple& strides() const;
	const CoordinateOrder& coordinateOrder() const;
	const bool& isSimple() const; 
	template<bool isConstLocal> bool overlaps(const View<T, isConstLocal>&);

	// element access
	template<class U> reference operator()(U); 
	template<class U> reference operator()(U) const; 
	#ifndef HAVE_CPP0X_VARIADIC_TEMPLATES
		reference operator()(const coordinate_type&,
			const coordinate_type&);
		reference operator()(const coordinate_type&,
			const coordinate_type&) const;
		reference operator()(const coordinate_type&,
			const coordinate_type&,
			const coordinate_type&);
		reference operator()(const coordinate_type&,
			const coordinate_type&,
			const coordinate_type&) const;
		reference operator()(const coordinate_type&,
			const coordinate_type&,
			const coordinate_type&,
			const coordinate_type&);
		reference operator()(const coordinate_type&,
			const coordinate_type&,
			const coordinate_type&,
			const coordinate_type&) const;
	#else
		template<typename... Args>
			reference operator()(const size_t &&, const Args && ...);
		reference operator()(const size_t &&);
		template<typename... Args>
			reference operator()(const size_t &&, const Args && ...) const;
		reference operator()(const size_t &&) const;
		private:
			template<typename... Args>
				size_t elementAccessHelper(const int &&, const size_t &&,
					const Args && ...);
			size_t elementAccessHelper(const int &&, const size_t &&);
			template<typename... Args>
				size_t elementAccessHelper(const int &&, const size_t &&, 
					const Args && ...) const;
			size_t elementAccessHelper(const int &&, const size_t &&) const;
		public:
	#endif

	// sub-views
	template<class BaseIterator, class ShapeIterator>
		void view(BaseIterator, ShapeIterator, View<T, isConst>&) const;
	template<class BaseIterator, class ShapeIterator>
		void view(BaseIterator, ShapeIterator, const CoordinateOrder&,
			View<T, isConst>&) const;
    template<class BaseIterator, class ShapeIterator>
        View<T, isConst> view(BaseIterator, ShapeIterator) const;
    template<class BaseIterator, class ShapeIterator>
        View<T, isConst> view(BaseIterator, ShapeIterator,
			const CoordinateOrder&) const;
	template<class BaseIterator, class ShapeIterator>
		void constView(BaseIterator, ShapeIterator, View<T, true>&) const;
	template<class BaseIterator, class ShapeIterator>
		void constView(BaseIterator, ShapeIterator, const CoordinateOrder&,
			View<T, true>&) const;
    template<class BaseIterator, class ShapeIterator>
        View<T, true> constView(BaseIterator, ShapeIterator) const;
    template<class BaseIterator, class ShapeIterator>
        View<T, true> constView(BaseIterator, ShapeIterator, 
			const CoordinateOrder&) const;
	#ifdef HAVE_CPP0X_INITIALIZER_LISTS
		void view(std::initializer_list<size_t>,
			std::initializer_list<size_t>, View<T, isConst>&) const;
		void view(std::initializer_list<size_t>,
			std::initializer_list<size_t>, const CoordinateOrder&,
			View<T, isConst>&) const;
		void constView(std::initializer_list<size_t>,
			std::initializer_list<size_t>, View<T, true>&) const;
		void constView(std::initializer_list<size_t>,
			std::initializer_list<size_t>, const CoordinateOrder&,
			View<T, true>&) const;
	#endif
	#ifdef MARRAY_COMPATIBILITY
		template<class BaseIterator, class EndIterator>
			View<T, isConst> subarray(BaseIterator, EndIterator) const;
	#endif

	// iterator access
	iterator begin();
	iterator end();
	const_iterator begin() const;
	const_iterator end() const;
	reverse_iterator rbegin();
	reverse_iterator rend();
	const_reverse_iterator rbegin() const;
	const_reverse_iterator rend() const;

	// coordinate transformation
	template<class ShapeIterator>
		void reshape(ShapeIterator, ShapeIterator);
	template<class CoordinateIterator>
		void permute(CoordinateIterator);
	void transpose(const coordinate_type&, const coordinate_type&);
	void transpose();
	void shift(int);
	void squeeze();

	template<class ShapeIterator>
		View<T, isConst> reshapedView(ShapeIterator, ShapeIterator) const;
	template<class CoordinateIterator>
		View<T, isConst> permutedView(CoordinateIterator) const;
	View<T, isConst> transposedView(const coordinate_type&, const coordinate_type&) const;
	View<T, isConst> transposedView() const;
	View<T, isConst> shiftedView(int) const;
	View<T, isConst> boundView(const dimension_type&, const coordinate_type& = 0) const;
	View<T, isConst> squeezedView() const;

	#ifdef HAVE_CPP0X_INITIALIZER_LISTS
		void reshape(std::initializer_list<size_t>);
		void permute(std::initializer_list<size_t>);

		View<T, isConst> reshapedView(std::initializer_list<size_t>) const;
		View<T, isConst> permutedView(std::initializer_list<size_t>) const;
	#endif

	#ifdef MARRAY_COMPATIBILITY
		template<class CoordinateIterator>
			View<T, isConst> permuteDimensions(CoordinateIterator) const;
		View<T, isConst> shiftDimensions(int) const;
	#endif

	// coordinates, index, offset conversion
	template<class CoordinateIterator>
		void coordinatesToIndex(CoordinateIterator, size_t&) const;
	template<class CoordinateIterator>
		void coordinatesToOffset(CoordinateIterator, size_t&) const;
	template<class CoordinateIterator>
		void indexToCoordinates(size_t, CoordinateIterator) const;
	void indexToOffset(const size_t&, size_t&) const;
	#ifdef HAVE_CPP0X_INITIALIZER_LISTS
		void coordinatesToIndex(std::initializer_list<size_t>,
			size_t&) const;
		void coordinatesToOffset(std::initializer_list<size_t>,
			size_t&) const;
	#endif

	// output as string
	std::string asString(const StringStyle& = MatrixStyle) const; 

private:
	void updateSimplicity();
	void testInvariant() const; 

	// shape and shape strides
	dimension_type dimension_;
	coordinate_tuple shape_;
	offset_tuple shapeStrides_;
		// Shape strides. Intended redundancy: shapeStrides_ could be
		// computed from shape_ and coordinateOrder_
	size_t size_;
		// size_ is the maximum index minus one, which may be unequal
		// to the length of data_ in memory, depending on the strides.
		// intended redundancy: size_ could be computed from shape_
	CoordinateOrder coordinateOrder_;

	// data and memory address
	pointer data_;
	size_t offset_;
	offset_tuple strides_;
	bool isSimple_; 
		// simple array: an array which is unstrided (i.e. the strides
		// equal the shape strides) and has a zero offset (cf. the
		// function testInvariant of View for the formal definition).

friend class View<T, false>;
friend class View<T, true>;
friend struct AssignmentOperatorHelper<T, true>;
friend struct AssignmentOperatorHelper<T, false>;
friend struct AccessOperatorHelper<true>;
friend struct AccessOperatorHelper<false>;
friend class Marray<T>;
friend class Vector<T>;
friend class Matrix<T>;
};

// arithmetic operators

template<class T>
	View<T, false>& operator+=(View<T, false>&, const T&);
template<class T1, class T2, bool isConst>
	View<T1, false>& operator+=(View<T1, false>&, const View<T2, isConst>&);
template<class T>
	View<T, false>& operator++(View<T, false>&); // prefix
template<class T>
	Marray<T> operator++(Marray<T>&, int); // postfix
template<class T, bool isConst>
	Marray<T> operator+(const View<T, isConst>&, const T&);
template<class T, bool isConst>
	Marray<T> operator+(const T&, const View<T, isConst>&);
template<class T, bool isConst1, bool isConst2>
	Marray<T> operator+(const View<T, isConst1>&, const View<T, isConst2>&);
template<class T, bool isConst>
	Marray<T> operator+(const View<T, isConst>&); // unary

template<class T>
	View<T, false>& operator-=(View<T, false>&, const T&);
template<class T1, class T2, bool isConst>
	View<T1, false>& operator-=(View<T1, false>&, const View<T2, isConst>&);
template<class T>
	View<T, false>& operator--(View<T, false>&); // prefix
template<class T>
	Marray<T> operator--(Marray<T>&, int); // postfix
template<class T, bool isConst>
	Marray<T> operator-(const View<T, isConst>&, const T&);
template<class T, bool isConst>
	Marray<T> operator-(const T&, const View<T, isConst>&);
template<class T, bool isConst1, bool isConst2>
	Marray<T> operator-(const View<T, isConst1>&, const View<T, isConst2>&);
template<class T, bool isConst>
	Marray<T> operator-(const View<T, isConst>&); // unary

template<class T>
	View<T, false>& operator*=(View<T, false>&, const T&);
template<class T1, class T2, bool isConst>
	View<T1, false>& operator*=(View<T1, false>&, const View<T2, isConst>&);
template<class T, bool isConst>
	Marray<T> operator*(const View<T, isConst>&, const T&);
template<class T, bool isConst>
	Marray<T> operator*(const T&, const View<T, isConst>&);
template<class T, bool isConst1, bool isConst2>
	Marray<T> operator*(const View<T, isConst1>&, const View<T, isConst2>&);

template<class T>
	View<T, false>& operator/=(View<T, false>&, const T&);
template<class T1, class T2, bool isConst>
	View<T1, false>& operator/=(View<T1, false>&, const View<T2, isConst>&);
template<class T, bool isConst>
	Marray<T> operator/(const View<T, isConst>&, const T&);
template<class T, bool isConst>
	Marray<T> operator/(const T&, const View<T, isConst>&);
template<class T, bool isConst1, bool isConst2>
	Marray<T> operator/(const View<T, isConst1>&, const View<T, isConst2>&);

/// STL compliant random access iterator for View and Marray.
/// 
/// In addition to the STL iterator interface, the member functions
/// hasMore(), index(), and coordinate() are defined.
template<class T, bool isConst>
class Iterator
{
public:
	// STL random access iterator typedefs
	typedef typename std::random_access_iterator_tag iterator_category;
	typedef T value_type;
	typedef ptrdiff_t difference_type;
	typedef typename IfBool<isConst, const T*, T*>::type pointer;
	typedef typename IfBool<isConst, const T&, T&>::type reference;

	// non-standard typedefs
	typedef size_t coordinate_type;
	typedef	typename IfBool<isConst, const View<T, true>*,
		View<T, false>*>::type view_pointer;
	typedef	typename IfBool<isConst, const View<T, true>&,
		View<T, false>&>::type view_reference;

	// construction
	Iterator();
	Iterator(const View<T, false>&, const size_t& = 0); 
	Iterator(View<T, false>&, const size_t& = 0); 
	Iterator(const View<T, true>&, const size_t& = 0);
	Iterator(const Iterator<T, false>&);
		// conversion from mutable to const

	// STL random access iterator operations
	reference operator*() const;
	pointer operator->() const;
	reference operator[](const size_t&) const;
	Iterator<T, isConst>& operator+=(const difference_type&);
	Iterator<T, isConst>& operator-=(const difference_type&);
	Iterator<T, isConst>& operator++(); // prefix
	Iterator<T, isConst>& operator--(); // prefix
	Iterator<T, isConst> operator++(int); // postfix
	Iterator<T, isConst> operator--(int); // postfix
	Iterator<T, isConst> operator+(const difference_type&) const;
	Iterator<T, isConst> operator-(const difference_type&) const;
	template<bool isConstLocal>
		difference_type operator-(const Iterator<T, isConstLocal>&) const;
	template<bool isConstLocal>
		bool operator==(const Iterator<T, isConstLocal>&) const;
	template<bool isConstLocal>
		bool operator!=(const Iterator<T, isConstLocal>&) const;
	template<bool isConstLocal>
		bool operator<(const Iterator<T, isConstLocal>&) const;
	template<bool isConstLocal>
		bool operator>(const Iterator<T, isConstLocal>&) const;
	template<bool isConstLocal>
		bool operator<=(const Iterator<T, isConstLocal>&) const;
	template<bool isConstLocal>
		bool operator>=(const Iterator<T, isConstLocal>&) const;

	// operations beyond the STL standard
	bool hasMore() const; 
	size_t index() const;
	template<class CoordinateIterator>
		void coordinate(CoordinateIterator) const;

private:
	// attributes
	view_pointer view_; 
	size_t index_;

friend class Marray<T>;
friend class Iterator<T, !isConst>; // for comparison operators
};

/// Runtime-Flexible multi-dimensional array.
template<class T> 
class Marray
: public View<T, false>
{
public:
	typedef View<T, false> base;
	typedef typename base::value_type value_type;
	typedef typename base::pointer pointer;
	typedef typename base::const_pointer const_pointer;
	typedef typename base::reference reference;
	typedef typename base::const_reference const_reference;
	typedef typename base::dimension_type dimension_type;
	typedef typename base::coordinate_type coordinate_type;
	typedef typename base::iterator iterator;
	typedef typename base::reverse_iterator reverse_iterator;
	typedef typename base::const_iterator const_iterator;
	typedef typename base::const_reverse_iterator const_reverse_iterator;
	typedef typename base::coordinate_tuple coordinate_tuple;
	typedef typename base::offset_tuple offset_tuple;
	#ifdef MARRAY_COMPATIBILITY
		typedef typename base::difference_type difference_type; ///< For compatibility with Vigra MultiArray
		typedef coordinate_type difference_type1; ///< For compatibility with Vigra MultiArray
	#endif

	// constructors and destructor
	Marray();
	Marray(const T&, const CoordinateOrder& = defaultOrder);
		// initialize as scalar
	Marray(const Marray<T>&);
	template<bool isConst>
		Marray(const View<T, isConst>&);
	template<class ShapeIterator>
		Marray(ShapeIterator, ShapeIterator, const T& = T(),
			const CoordinateOrder& = defaultOrder);
	template<class ShapeIterator>
		Marray(const InitializationSkipping&, ShapeIterator, ShapeIterator,
			const CoordinateOrder& = defaultOrder);
	~Marray();
	void assign();

	#ifdef HAVE_CPP0X_INITIALIZER_LISTS
		Marray(std::initializer_list<size_t>, const T& = T(),
			const CoordinateOrder& = defaultOrder);
	#endif

	// assignment operator
	Marray<T>& operator=(const Marray<T>&);
		// overwrite standard operator=
	template<bool isConst>
		Marray<T>& operator=(const View<T, isConst>&);

	// resize
	template<class ShapeIterator>
		void resize(ShapeIterator, ShapeIterator, const T& = T());
	#ifdef HAVE_CPP0X_INITIALIZER_LISTS
		void resize(std::initializer_list<size_t>, const T& = T());
	#endif

private:

	void testInvariant() const; 

friend class Vector<T>;
friend class Matrix<T>;
};

/// One-dimensional Marray.
template<class T> 
class Vector
: public Marray<T>
{
public:
	typedef Marray<T> base;
	typedef typename base::value_type value_type;
	typedef typename base::pointer pointer;
	typedef typename base::const_pointer const_pointer;
	typedef typename base::reference reference;
	typedef typename base::const_reference const_reference;
	typedef typename base::dimension_type dimension_type;
	typedef typename base::coordinate_type coordinate_type;
	typedef typename base::iterator iterator;
	typedef typename base::reverse_iterator reverse_iterator;
	typedef typename base::const_iterator const_iterator;
	typedef typename base::const_reverse_iterator const_reverse_iterator;
	typedef typename base::coordinate_tuple coordinate_tuple;
	typedef typename base::offset_tuple offset_tuple;
	#ifdef MARRAY_COMPATIBILITY
		typedef typename base::difference_type difference_type; ///< For compatibility with Vigra MultiArray
		typedef coordinate_type difference_type1; ///< For compatibility with Vigra MultiArray
	#endif

	// constructors and destructor
	Vector();
	template<bool isConst>
		Vector(const View<T, isConst>&);
	Vector(const size_t&, const T& = T());
	Vector(const InitializationSkipping&, const size_t&);
	#ifdef HAVE_CPP0X_INITIALIZER_LISTS
		Vector(std::initializer_list<T> list);
	#endif

	// assignment operator
	Vector<T>& operator=(const Vector<T>&);
		// overwrite standard operator=
	template<bool isConst>
		Vector<T>& operator=(const View<T, isConst>&);

	// element access
	T& operator[](const size_t&);
	const T& operator[](const size_t&) const;

	// reshape, resize
	void reshape(const size_t&);
	void resize(const size_t&, const T& = T());

private:
	void testInvariant() const;
};

/// Two-dimensional Marray.
template<class T> 
class Matrix
: public Marray<T>
{
public:
	typedef Marray<T> base;
	typedef typename base::value_type value_type;
	typedef typename base::pointer pointer;
	typedef typename base::const_pointer const_pointer;
	typedef typename base::reference reference;
	typedef typename base::const_reference const_reference;
	typedef typename base::dimension_type dimension_type;
	typedef typename base::coordinate_type coordinate_type;
	typedef typename base::iterator iterator;
	typedef typename base::reverse_iterator reverse_iterator;
	typedef typename base::const_iterator const_iterator;
	typedef typename base::const_reverse_iterator const_reverse_iterator;
	typedef typename base::coordinate_tuple coordinate_tuple;
	typedef typename base::offset_tuple offset_tuple;

	// constructors and destructor
	Matrix();
    template<bool isConst>
		Matrix(const View<T, isConst>&);
	Matrix(const size_t&, const size_t&, const T& = T(),
		const CoordinateOrder& = defaultOrder);
	Matrix(const InitializationSkipping&, const size_t&, const size_t&,
		const CoordinateOrder& = defaultOrder);

	// assignment operator
	Matrix<T>& operator=(const Matrix<T>&);
		// overwrite standard operator=
	template<bool isConst>
		Matrix<T>& operator=(const View<T, isConst>&);

	// resize and reshape
	void reshape(const size_t&, const size_t&);
	void resize(const size_t&, const size_t&, const T& = T());

private:
	void testInvariant() const;
};

// implementation of namespace functions

/// Compute a sequence of strides from a shape sequence.
template<class ShapeIterator, class StridesIterator>
void shapeStrides
(
	ShapeIterator begin,
	ShapeIterator end,
	StridesIterator strideBegin,
	const CoordinateOrder& coordinateOrder
) 
{
	Assert(NO_ARG_TEST || std::distance(begin, end) != 0);

	size_t dimension = std::distance(begin, end);
	ShapeIterator shapeIt;
	StridesIterator strideIt;
	if(coordinateOrder == FirstMajorOrder) {
		shapeIt = begin + (dimension-1);
		strideIt = strideBegin + (dimension-1);
	}
	else {
		shapeIt = begin;
		strideIt = strideBegin;
	}
	*strideIt = 1;
	for(size_t j=1; j<dimension; ++j) {
		size_t tmp = *strideIt;
		if(coordinateOrder == FirstMajorOrder) {
			--strideIt;
			(*strideIt) = tmp * (*shapeIt);
			--shapeIt;
		}
		else {
			++strideIt;
			(*strideIt) = tmp * (*shapeIt);
			++shapeIt;
		}
	}
}

// implementation of View

#ifdef HAVE_CPP0X_INITIALIZER_LISTS
/// Compute the index that corresponds to a sequence of coordinates.
///
/// \param coordinate Coordinate given as initializer list.
/// \param out Index (output)
/// \sa coordinatesToOffset(), indexToCoordinates(), and indexToOffset()
///
template<class T, bool isConst>
inline void
View<T, isConst>::coordinatesToIndex
(
	std::initializer_list<size_t> coordinate,
	size_t& out
) const 
{
	coordinatesToIndex(coordinate.begin(), out);
}
#endif

/// Compute the index that corresponds to a sequence of coordinates.
///
/// \param it An iterator to the beginning of the coordinate sequence.
/// \param out Index (output)
/// \sa coordinatesToOffset(), indexToCoordinates(), and indexToOffset()
///
template<class T, bool isConst> 
template<class CoordinateIterator>
inline void 
View<T, isConst>::coordinatesToIndex
(
	CoordinateIterator it,
	size_t& out
) const
{
	testInvariant();
	out = 0;
	for(size_t j=0; j<this->dimension_; ++j) {
		Assert(NO_ARG_TEST || coordinate_type(*it) < shape_[j]);
		out += static_cast<size_t>(*it) * shapeStrides_[j];
		++it;
	}
}

#ifdef HAVE_CPP0X_INITIALIZER_LISTS
/// Compute the offset that corresponds to a sequence of coordinates.
///
/// \param it An iterator to the beginning of the coordinate sequence.
/// \param out Index (output)
/// \sa coordinatesToIndex(), indexToCoordinates(), and indexToOffset()
///
template<class T, bool isConst> 
inline void
View<T, isConst>::coordinatesToOffset
(
	std::initializer_list<size_t> coordinate,
	size_t& out
) const
{
	coordinatesToOffset(coordinate.begin(), out);
}
#endif

/// Compute the offset that corresponds to a sequence of coordinates.
///
/// \param it An iterator to the beginning of the coordinate sequence.
/// \param out Index (output)
/// \sa coordinatesToIndex(), indexToCoordinates(), and indexToOffset()
///
template<class T, bool isConst> 
template<class CoordinateIterator>
inline void
View<T, isConst>::coordinatesToOffset
(
	CoordinateIterator it,
	size_t& out
) const
{
	testInvariant();
	out = offset_;
	for(size_t j=0; j<this->dimension_; ++j) {
		Assert(NO_ARG_TEST || coordinate_type(*it) < shape_[j]);
		out += static_cast<size_t>(*it) * strides_[j];
		++it;
	}
}

/// Compute the coordinate sequence that corresponds to an index.
///
/// \param index Index
/// \param outit An iterator into a container into which the coordinate
/// sequence is to be written (output).
/// \sa coordinatesToIndex(), coordinatesToOffset(), and indexToOffset()
///
template<class T, bool isConst>
template<class CoordinateIterator>
inline void
View<T, isConst>::indexToCoordinates
(
	size_t index, // copy to work on
	CoordinateIterator outit
) const
{
	testInvariant();
	Assert(NO_DEBUG || this->dimension_ > 0);
	Assert(NO_ARG_TEST || index < this->size_);

	if(coordinateOrder_ == FirstMajorOrder) {
		for(size_t j=0; j<this->dimension_; ++j) {
			*outit = coordinate_type(index / shapeStrides_[j]);
			index = index % shapeStrides_[j];
			++outit;
		}
	}
	else { // last major order
		size_t j = this->dimension_-1;
		outit += j;
		for(;;) {
			*outit = coordinate_type(index / shapeStrides_[j]);
			index = index % shapeStrides_[j];
			if(j == 0) {
				break;
			}
			else {
				--outit;
				--j;
			}
		}
	}
}

/// Compute the offset that corresponds to an index.
///
/// \param index Index.
/// \param out Offset (output).
/// \sa coordinatesToIndex(), coordinatesToOffset(), and indexToCoordinates()
///
template<class T, bool isConst> 
inline void
View<T, isConst>::indexToOffset
(
	const size_t& index,
	size_t& out
) const
{
	testInvariant();
	Assert(NO_ARG_TEST || index < this->size_); 
	
	if(isSimple_) {
		out = index;
	}
	else {
		coordinate_tuple coordinates(shape_.size()); 
		indexToCoordinates(index, coordinates.begin());
		coordinatesToOffset(coordinates.begin(), out);
	}
}

/// Empty constructor.
///
/// The empty constructor sets the data pointer to 0.
/// It does not allocate memory for a scalar.
///
template<class T, bool isConst> 
inline
View<T, isConst>::View()
: dimension_(0),
  shape_(coordinate_tuple()),
  shapeStrides_(offset_tuple()),
  size_(0),
  coordinateOrder_(defaultOrder),
  data_(0),
  offset_(0),
  strides_(offset_tuple()),
  isSimple_(true)
{
	testInvariant();
}

/// Construct a View to a scalar.
///
/// \param data Pointer to data.
///
template<class T, bool isConst> 
inline
View<T, isConst>::View
(
	pointer data
)
: dimension_(0),
  shape_(coordinate_tuple()),
  shapeStrides_(offset_tuple()),
  size_(1),
  coordinateOrder_(defaultOrder),
  data_(data),
  offset_(0),
  strides_(offset_tuple()),
  isSimple_(true)
{
	testInvariant();
}

/// Cast from a View on mutable data to a View on constant data.
///
template<class T, bool isConst> 
inline
View<T, isConst>::View
(
	const View<T, false>& in
)
: dimension_(in.dimension_),
  shape_(in.shape_),
  shapeStrides_(in.shapeStrides_),
  size_(in.size_),
  coordinateOrder_(in.coordinateOrder_),
  data_(in.data_),
  offset_(in.offset_),
  strides_(in.strides_),
  isSimple_(in.isSimple_)
{
	testInvariant();
}

#ifdef HAVE_CPP0X_INITIALIZER_LISTS
/// Construct unstrided View with zero offset
/// 
/// \param shape Shape initializer list.
/// \param data Pointer to data.
/// \param externalCoordinateOrder Flag specifying the order
/// of coordinates based on which the strides are computed.
/// \param internalCoordinateOrder Flag specifying the order
/// of coordinates used for scalar indexing and iterators.
///

template<class T, bool isConst>
inline
View<T, isConst>::View
(
	std::initializer_list<size_t> shape,
	pointer data,
	const CoordinateOrder& externalCoordinateOrder,
	const CoordinateOrder& internalCoordinateOrder
)
{
	assign(shape.begin(), shape.end(), data, externalCoordinateOrder,
		internalCoordinateOrder);
}
#endif

/// Construct unstrided View with zero offset
/// 
/// \param begin Iterator to the beginning of a sequence that
/// defines the shape.
/// \param end Iterator to the end of this sequence.
/// \param data Pointer to data.
/// \param externalCoordinateOrder Flag specifying the order
/// of coordinates based on which the strides are computed.
/// \param internalCoordinateOrder Flag specifying the order
/// of coordinates used for scalar indexing and iterators.
///
template<class T, bool isConst> 
template<class ShapeIterator>
inline
View<T, isConst>::View
(
	ShapeIterator begin,
	ShapeIterator end,
	pointer data,
	const CoordinateOrder& externalCoordinateOrder,
	const CoordinateOrder& internalCoordinateOrder
)
{
	assign(begin, end, data, externalCoordinateOrder,
		internalCoordinateOrder);
}

#ifdef HAVE_CPP0X_INITIALIZER_LISTS
/// Construct strided View, with or without offset
/// 
/// \param shape Shape initializer list.
/// \param strides Strides initializer list.
/// \param data Pointer to data.
/// \param offset Offset.
/// \param internalCoordinateOrder Flag specifying the order
/// of coordinates used for scalar indexing and iterators.
///
template<class T, bool isConst> 
inline
View<T, isConst>::View
(
	std::initializer_list<size_t> shape,
	std::initializer_list<size_t> strides,
	pointer data,
	const size_t& offset,
	const CoordinateOrder& internalCoordinateOrder
)
{
	assign(shape.begin(), shape.end(), strides.begin(), data, offset,
		internalCoordinateOrder);
}
#endif

/// Construct strided View, with or without offset
/// 
/// \param begin Iterator to the beginning of a sequence that
/// defines the shape.
/// \param end Iterator to the end of this sequence.
/// \param it Iterator to the beginning of a sequence that
/// defines the strides.
/// \param data Pointer to data.
/// \param offset Offset.
/// \param internalCoordinateOrder Flag specifying the order
/// of coordinates used for scalar indexing and iterators.
///
template<class T, bool isConst> 
template<class ShapeIterator, class StrideIterator>
inline
View<T, isConst>::View
(
	ShapeIterator begin,
	ShapeIterator end,
	StrideIterator it,
	pointer data,
	const size_t& offset,
	const CoordinateOrder& internalCoordinateOrder
)
{
	assign(begin, end, it, data, offset, internalCoordinateOrder);
}

#ifdef HAVE_CPP0X_INITIALIZER_LISTS
/// Initialize unstrided View with zero offset
/// 
/// \param shape Shape initializer list.
/// \param data Pointer to data.
/// \param externalCoordinateOrder Flag specifying the order
/// of coordinates based on which the strides are computed.
/// \param internalCoordinateOrder Flag specifying the order
/// of coordinates used for scalar indexing and iterators.
///
template<class T, bool isConst> 
void
View<T, isConst>::assign
(
	std::initializer_list<size_t> shape,
	pointer data,
	const CoordinateOrder& externalCoordinateOrder,
	const CoordinateOrder& internalCoordinateOrder
)
{
	assign(shape.begin(), shape.end(), data, externalCoordinateOrder,
		internalCoordinateOrder);
}
#endif

/// Clear View.
///
/// Leaves the View in the same state as if the empty constructor
/// had been called.
///
/// \sa View()
///
template<class T, bool isConst> 
inline void
View<T, isConst>::assign()
{
  this->dimension_ = 0;
  this->shape_ = coordinate_tuple();
  this->shapeStrides_ = offset_tuple();
  this->size_ = 0;
  this->coordinateOrder_ = defaultOrder;
  this->offset_ = 0;
  this->strides_ = offset_tuple();
  this->isSimple_ = true;
  this->data_ = 0;
	testInvariant();
}

/// Initialize unstrided View with zero offset
/// 
/// \param begin Iterator to the beginning of a sequence that
/// defines the shape.
/// \param end Iterator to the end of this sequence.
/// \param data Pointer to data.
/// \param externalCoordinateOrder Flag specifying the order
/// of coordinates based on which the strides are computed.
/// \param internalCoordinateOrder Flag specifying the order
/// of coordinates used for scalar indexing and iterators.
///
template<class T, bool isConst> 
template<class ShapeIterator>
void
View<T, isConst>::assign
(
	ShapeIterator begin,
	ShapeIterator end,
	pointer data,
	const CoordinateOrder& externalCoordinateOrder,
	const CoordinateOrder& internalCoordinateOrder
)
{
	// invariant is not tested as a pre-condition
	// because assign is called from constructors

	dimension_ = dimension_type(end - begin);
	shape_ = coordinate_tuple(dimension_);
	shapeStrides_ = offset_tuple(dimension_);
	data_ = data;
	offset_ = 0;
	strides_ = offset_tuple(dimension_);
	coordinateOrder_ = internalCoordinateOrder;
	size_ = 1;

	if(dimension_ == 0) { // if the array is a scalar
		isSimple_ = true;
	}
	else { // if the array is not a scalar
		isSimple_ = (externalCoordinateOrder == internalCoordinateOrder);
		for(size_t j=0; j<dimension_; ++j) {
			shape_[j] = coordinate_type(*begin);
			size_ *= (size_t)(*begin);
			++begin;
		}
		shapeStrides(shape_.begin(), shape_.end(), strides_.begin(),
			externalCoordinateOrder);
		shapeStrides(shape_.begin(), shape_.end(), shapeStrides_.begin(),
			internalCoordinateOrder);
	}

	testInvariant();
}

#ifdef HAVE_CPP0X_INITIALIZER_LISTS
/// Initialize strided View, with or without offset
/// 
/// \param shape Shape initializer list.
/// \param strides Strides initialier list.
/// defines the strides.
/// \param data Pointer to data.
/// \param offset Offset.
/// \param internalCoordinateOrder Flag specifying the order
/// of coordinates used for scalar indexing and iterators.
///
template<class T, bool isConst> 
void
View<T, isConst>::assign
(
	std::initializer_list<size_t> shape,
	std::initializer_list<size_t> strides,
	pointer data,
	const size_t& offset,
	const CoordinateOrder& internalCoordinateOrder
)
{
	assign(shape.begin(), shape.end(), strides.begin(), data, offset,
		internalCoordinateOrder);
}
#endif

/// Initialize strided View, with or without offset
/// 
/// \param begin Iterator to the beginning of a sequence that
/// defines the shape.
/// \param end Iterator to the end of this sequence.
/// \param it Iterator to the beginning of a sequence that
/// defines the strides.
/// \param data Pointer to data.
/// \param offset Offset.
/// \param internalCoordinateOrder Flag specifying the order
/// of coordinates used for scalar indexing and iterators.
///
template<class T, bool isConst> 
template<class ShapeIterator, class StrideIterator>
void
View<T, isConst>::assign
(
	ShapeIterator begin,
	ShapeIterator end,
	StrideIterator it,
	pointer data,
	const size_t& offset,
	const CoordinateOrder& internalCoordinateOrder
)
{
	// invariant is not tested as a pre-condition
	// because assign is called from constructors

	dimension_ = std::distance(begin, end);
	shape_ = coordinate_tuple(dimension_);
	shapeStrides_ = offset_tuple(dimension_);
	data_ = data;
	offset_ = offset;
	strides_ = offset_tuple(dimension_);
	coordinateOrder_ = internalCoordinateOrder;
	size_ = 1;

	if(dimension_ == 0) {
		if(offset_ == 0) {
			isSimple_ = true;
		}
		else {
			isSimple_ = false;
		}
	}
	else {
		for(size_t j=0; j<dimension_; ++j) {
			shape_[j] = *begin;
			size_ *= *begin;
			++begin;
			strides_[j] = *it;
			++it;
		}
		shapeStrides(shape_.begin(), shape_.end(), shapeStrides_.begin(),
			internalCoordinateOrder);
		updateSimplicity();
	}

	testInvariant();
}

// \cond suppress doxygen
template<>
struct AccessOperatorHelper<true>
{
	// access by scalar index
	template<class T, class U, bool isConst>
	static typename View<T, isConst>::reference
	execute(const View<T, isConst>& v, const U& index)
	{
		v.testInvariant();
		Assert(NO_DEBUG || v.data_ != 0);
		if(v.dimension_ == 0) {
			Assert(NO_ARG_TEST || index == 0);
			return v.data_[v.offset_];
		}
		else {
			size_t offset;
			v.indexToOffset(index, offset);
			return v.data_[offset];
		}
	}
};

template<>
struct AccessOperatorHelper<false>
{
	// access by iterator
	template<class T, class U, bool isConst>
	static typename View<T, isConst>::reference
	execute(const View<T, isConst>& v, U it)
	{
		v.testInvariant();
		Assert(NO_DEBUG || v.data_ != 0);
		size_t offset;
		v.coordinatesToOffset(it, offset);
		return v.data_[offset];
	}
};
// \endcond end suppress doxygen

/// Reference data.
///
/// \param u If u is an integer type, scalar indexing is performed.
/// Otherwise, it is assumed that u is an iterator to the beginning
/// of a coordinate sequence. 
/// \return Reference to the entry at u.
///
template<class T, bool isConst> 
template<class U>
inline typename View<T, isConst>::reference
View<T, isConst>::operator()
(
	U u
)
{
	return AccessOperatorHelper<std::numeric_limits<U>::is_integer>::execute(*this, u);
}

/// Reference data.
///
/// \param u If u is an integer type, scalar indexing is performed.
/// Otherwise, it is assumed that u is an iterator to the beginning
/// of a coordinate sequence. 
/// \return Reference to the entry at u.
///
template<class T, bool isConst> 
template<class U>
inline typename View<T, isConst>::reference
View<T, isConst>::operator()
(
	U u
) const
{
	return AccessOperatorHelper<std::numeric_limits<U>::is_integer>::execute(*this, u);
}

#ifndef HAVE_CPP0X_VARIADIC_TEMPLATES

/// Reference data in a 2-dimensional View by coordinates.
///
/// This function issues a runtime error if the View is not
/// 2-dimensional.
///
/// \param c0 1st coordinate. 
/// \param c1 2nd coordinate.
///
template<class T, bool isConst> 
inline typename View<T, isConst>::reference
View<T, isConst>::operator()
(
	const coordinate_type& c0,
    const coordinate_type& c1
)
{
	testInvariant();
	Assert(NO_DEBUG || (data_ != 0 && dimension_ == 2));
	Assert(NO_ARG_TEST || (c0 < shape_[0] && c1 < shape_[1]));

	return data_[offset_ + static_cast<size_t>(c0)*strides_[0]
		+ static_cast<size_t>(c1)*strides_[1] ];
}

/// Reference data in a 2-dimensional View by coordinates.
///
/// This function issues a runtime error if the View is not
/// 2-dimensional.
///
/// \param c0 1st coordinate. 
/// \param c1 2nd coordinate.
///
template<class T, bool isConst> 
inline typename View<T, isConst>::reference
View<T, isConst>::operator()
(
	const coordinate_type& c0,
	const coordinate_type& c1
) const
{
	testInvariant();
	Assert(NO_DEBUG || (data_ != 0 && dimension_ == 2));
	Assert(NO_ARG_TEST || (c0 < shape_[0] && c1 < shape_[1]));

	return data_[ offset_ + static_cast<size_t>(c0)*strides_[0]
		+ static_cast<size_t>(c1)*strides_[1] ];
}

/// Reference data in a 3-dimensional View by coordinates.
///
/// This function issues a runtime error if the View is not
/// 3-dimensional.
///
/// \param c0 1st coordinate. 
/// \param c1 2nd coordinate.
/// \param c2 3rd coordinate.
///
template<class T, bool isConst> 
inline typename View<T, isConst>::reference
View<T, isConst>::operator()
(
	const coordinate_type& c0,
	const coordinate_type& c1,
	const coordinate_type& c2
)
{
	testInvariant();
	Assert(NO_DEBUG || (data_ != 0 && dimension_ == 3));
	Assert(NO_ARG_TEST || (c0 < shape_[0] && c1 < shape_[1]
		&& c2 < shape_[2]));

	return data_[offset_ + static_cast<size_t>(c0)*strides_[0]
		+ static_cast<size_t>(c1)*strides_[1] + static_cast<size_t>(c2)*strides_[2] ];
}

/// Reference data in a 3-dimensional View by coordinates.
///
/// This function issues a runtime error if the View is not
/// 3-dimensional.
///
/// \param c0 1st coordinate. 
/// \param c1 2nd coordinate.
/// \param c2 3rd coordinate.
///
template<class T, bool isConst> 
inline typename View<T, isConst>::reference
View<T, isConst>::operator()
(
	const coordinate_type& c0,
	const coordinate_type& c1,
	const coordinate_type& c2
) const
{
	testInvariant();
	Assert(NO_DEBUG || (data_ != 0 && dimension_ == 3));
	Assert(NO_ARG_TEST || (c0 < shape_[0] && c1 < shape_[1]
		&& c2 < shape_[2]));

	return data_[ offset_ + static_cast<size_t>(c0)*strides_[0]
		+ static_cast<size_t>(c1)*strides_[1] + static_cast<size_t>(c2)*strides_[2] ];
}

/// Reference data in a 4-dimensional View by coordinates.
///
/// This function issues a runtime error if the View is not
/// 4-dimensional.
///
/// \param c0 1st coordinate. 
/// \param c1 2nd coordinate.
/// \param c2 3rd coordinate.
/// \param c3 4th coordinate.
///
template<class T, bool isConst> 
inline typename View<T, isConst>::reference
View<T, isConst>::operator()
(
	const coordinate_type& c0,
	const coordinate_type& c1,
	const coordinate_type& c2,
	const coordinate_type& c3
)
{
	testInvariant();
	Assert(NO_DEBUG || (data_ != 0 && dimension_ == 4));
	Assert(NO_ARG_TEST || (c0 < shape_[0] && c1 < shape_[1]
		&& c2 < shape_[2] && c3 < shape_[3]));

	return data_[offset_ + static_cast<size_t>(c0)*strides_[0]
		+ static_cast<size_t>(c1)*strides_[1] + static_cast<size_t>(c2)*strides_[2]
		+ static_cast<size_t>(c3)*strides_[3] ];
}

/// Reference data in a 4-dimensional View by coordinates.
///
/// This function issues a runtime error if the View is not
/// 4-dimensional.
///
/// \param c0 1st coordinate. 
/// \param c1 2nd coordinate.
/// \param c2 3rd coordinate.
/// \param c3 4th coordinate.
///
template<class T, bool isConst> 
inline typename View<T, isConst>::reference
View<T, isConst>::operator()
(
	const coordinate_type& c0,
	const coordinate_type& c1,
	const coordinate_type& c2,
	const coordinate_type& c3
) const
{
	testInvariant();
	Assert(NO_DEBUG || (data_ != 0 && dimension_ == 4));
	Assert(NO_ARG_TEST || (c0 < shape_[0] && c1 < shape_[1]
		&& c2 < shape_[2] && c3 < shape_[3]));

	return data_[ offset_ + static_cast<size_t>(c0)*strides_[0]
		+ static_cast<size_t>(c1)*strides_[1] + static_cast<size_t>(c2)*strides_[2]
		+ static_cast<size_t>(c3)*strides_[3] ];
}

#else

/// Reference data by a sequence of coordinates (C++0x)
///
/// This function takes an arbitrary number of coordinates. If this
/// number does not match with the dimension, a runtime error is issued.
///
template<class T, bool isConst>
template<typename... Args>
inline typename View<T, isConst>::reference
View<T, isConst>::operator()(const size_t && value, const Args && ... args)
{
	testInvariant();
	
	//one argument unpacked, so the size of args is one less than the dimension
	Assert( NO_DEBUG || ( data_ != 0 && sizeof...(args)+1 == dimension() ) );
	
	//We pass the rest of the arguments to the elementAccessHelper helper function,
	//which calculates the rest of the index expression used to access the
	//linear data_ array at the correct position
	return data_[ offset_ + strides_[0]*value + elementAccessHelper(sizeof...(args)+1, args...) ];
}

/// Reference data by a sequence of coordinates (C++0x)
///
/// This function takes an arbitrary number of coordinates. If this
/// number does not match with the dimension, a runtime error is issued.
///
template<class T, bool isConst>
inline typename View<T, isConst>::reference
View<T, isConst>::operator()(const size_t && value)
{
	testInvariant();
	
	if(dimension_ == 0) {
		Assert(NO_ARG_TEST || value == 0);
		return data_[offset_];
	}
	else {
		size_t offset;
		indexToOffset(value, offset);
		return data_[offset];
	}
}

// the element access helper recursively computes the remaining
// expression for the correct array offset.
//
// We pass in the dimension to allow the compiler to completely
// reduce the expression to the above non-c++0x definitions. The
// dimension is known at compile time through the number of
// arguments passed. We therefore should not use the dimension()
// call.
//
template<class T, bool isConst>
template<typename... Args>
inline size_t
View<T, isConst>::elementAccessHelper(const int && Dim, const size_t && value, const Args && ... args)
{
	Assert(NO_ARG_TEST || (value < shape_[Dim-1-sizeof...(args)] ) );
		
	return value*strides_[Dim-1-sizeof...(args)] + elementAccessHelper(Dim, args...); 
}

// basis case of recursion: only one argument left
template<class T, bool isConst>
inline size_t
View<T, isConst>::elementAccessHelper(const int && Dim, const size_t && value)
{
	Assert(NO_ARG_TEST || (value < shape_[Dim-1] ) );
	
	return strides_[Dim-1]*value;
}

template<class T, bool isConst>
template<typename... Args>
inline typename View<T, isConst>::reference
View<T, isConst>::operator()(const size_t && value, const Args && ... args) const
{
	testInvariant();
	
	Assert( NO_DEBUG || ( data_ != 0 && sizeof...(args)+1 == dimension() ) );
	
	return data_[ offset_ + strides_[0]*static_cast<size_t>(value) 
		+ static_cast<size_t>(elementAccessHelper(sizeof...(args)+1, args...)) ];
}

template<class T, bool isConst>
inline typename View<T, isConst>::reference
View<T, isConst>::operator()(const size_t && value) const
{
	testInvariant();
	
	if(dimension_ == 0) {
		Assert(NO_ARG_TEST || value == 0);
		return data_[offset_];
	}
	else {
		size_t offset;
		indexToOffset(value, offset);
		return data_[offset];
	}
}

template<class T, bool isConst>
template<typename... Args>
inline size_t
View<T, isConst>::elementAccessHelper(const int && Dim, const size_t && value, const Args && ... args) const
{
	Assert(NO_ARG_TEST || (value < shape_[Dim-1-sizeof...(args)] ) );
	
	return value*strides_[Dim-1-sizeof...(args)] + elementAccessHelper(Dim, args...); 
}

template<class T, bool isConst>
inline size_t
View<T, isConst>::elementAccessHelper(const int && Dim, const size_t && value) const
{
	Assert(NO_ARG_TEST || (value < shape_[Dim-1] ) );
	
	return strides_[Dim-1]*value;
}

#endif // #ifndef HAVE_CPP0X_VARIADIC_TEMPLATES

/// Get the number of data items.
///
/// \return Size.
///
template<class T, bool isConst> 
inline const size_t&
View<T, isConst>::size() const
{
	return size_;
}

/// Get the dimension.
///
/// \return Dimension.
///
template<class T, bool isConst> 
inline const typename View<T, isConst>::dimension_type&
View<T, isConst>::dimension() const
{
  Assert(NO_DEBUG || this->data_ != 0); 
    // the dimension is not well-defined if this->data_ == 0
	return dimension_;
}

/// Get the shape in one dimension.
///
/// \param dimension. Dimension
/// \return Shape in that dimension.
///
template<class T, bool isConst> 
inline const typename View<T, isConst>::coordinate_type&
View<T, isConst>::shape
(
	const typename View<T, isConst>::dimension_type& dimension
) const
{
	testInvariant();
	Assert(NO_DEBUG || data_ != 0);
	Assert(NO_ARG_TEST || dimension < dimension_);
	return shape_[dimension];
}

/// Get a reference to the shape.
///
/// \return Reference to the shape coordinate tuple.
///
template<class T, bool isConst> 
inline const typename View<T, isConst>::coordinate_tuple&
View<T, isConst>::shape() const
{
	return shape_;
}

/// Get the strides in one dimension.
///
/// \param dimension. Dimension
/// \return Stride in that dimension.
///
template<class T, bool isConst> 
inline const size_t&
View<T, isConst>::strides
(
	const typename View<T, isConst>::dimension_type& dimension
) const
{
	testInvariant();
	Assert(NO_DEBUG || data_ != 0);
	Assert(NO_ARG_TEST || dimension < dimension_);
	return strides_[dimension];
}

/// Get a reference to the strides.
///
/// \return Reference to the strides offset tuple.
///
template<class T, bool isConst> 
inline const typename View<T, isConst>::offset_tuple&
View<T, isConst>::strides() const
{
	return strides_;
}

/// Get the coordinate order used for scalar indexing and iterators.
///
/// \return CoordinateOrder. enum: FirstMajorOrder, LastMajorOrder
///
template<class T, bool isConst> 
inline const CoordinateOrder&
View<T, isConst>::coordinateOrder() const
{
	testInvariant();
	return coordinateOrder_;
}

/// Determine whether the Marray is unstrided and has zero offset.
///
/// \return bool.
///
template<class T, bool isConst> 
inline const bool& 
View<T, isConst>::isSimple() const
{
	testInvariant();
	return isSimple_;
}

// \cond suppress doxygen
/// AssignmentOperatorHelper class
///
/// This is a helper class for the class template View. It is used
/// to implement a different behavoir of the assignment operator of
/// View, depending on the template parameter isConst of View.
/// Cf. the definition of operator= in View for details.
///
template<class T> 
struct AssignmentOperatorHelper<T, false>
{
	/// from constant to mutable
	static void execute
	(
		const View<T, true>& from,
		View<T, false>& to
	)
	{
		if(!NO_ARG_TEST) {
			Assert(from.data_ != 0 && from.dimension_ == to.dimension_);
			for(size_t j=0; j<from.dimension_; ++j) {
				assert(from.shape(j) == to.shape(j));
			}
		}
		if(from.coordinateOrder_ == to.coordinateOrder_) {
			if(from.isSimple_ && to.isSimple_) {
				memcpy(to.data_, from.data_, (from.size_)*sizeof(T));
			}
			else {
				for(size_t j=0; j<from.size_; ++j) {
					to(j) = from(j);
				}
			}
		}
		else {
			std::vector<size_t> coords(from.dimension_);
			if(from.isSimple_) {
				for(size_t j=0; j<from.size_; ++j) {
					from.indexToCoordinates(j, coords.begin());
					to(coords.begin()) = from(j);
				}
			}
			else if(to.isSimple_) {
				for(size_t j=0; j<from.size_; ++j) {
					to.indexToCoordinates(j, coords.begin());
					to(j) = from(coords.begin());
				}
			}
			else {
				for(size_t j=0; j<from.size_; ++j) {
					to.indexToCoordinates(j, coords.begin());
					to(coords.begin()) = from(coords.begin());
				}
			}
		}
	}

	/// from mutable to mutable
	static void execute
	(
		const View<T, false>& from,
		View<T, false>& to
	)
	{
		if(&from != &to) { // no self-assignment
			if(to.data_ == 0) { // if the view 'to' is not initialized
				// initialize the view 'to' with source data
				to.dimension_ = from.dimension_;
				to.shape_ = from.shape_;
				to.shapeStrides_ = from.shapeStrides_;
				to.size_ = from.size_;
				to.data_ = from.data_; // copy pointer
				to.offset_ = from.offset_;
				to.strides_ = from.strides_;
				to.isSimple_ = from.isSimple_; 
				to.coordinateOrder_ = from.coordinateOrder_;
			}
			else { // if the view 'to' is initialized
				if(!NO_ARG_TEST) {
					Assert(from.data_ != 0 && from.dimension_ == to.dimension_);
					for(size_t j=0; j<from.dimension_; ++j) {
						assert(from.shape(j) == to.shape(j));
					}
				}
				if(from.coordinateOrder_ == to.coordinateOrder_) {
					if(from.isSimple_ && to.isSimple_) {
						memcpy(to.data_, from.data_, (from.size_)*sizeof(T));
					}
					else {
						for(size_t j=0; j<from.size_; ++j) {
							to(j) = from(j);
						}
					}
				}
				else {
					std::vector<size_t> coords(from.dimension_);
					if(from.isSimple_) {
						for(size_t j=0; j<from.size_; ++j) {
							from.indexToCoordinates(j, coords.begin());
							to(coords.begin()) = from(j);
						}
					}
					else if(to.isSimple_) {
						for(size_t j=0; j<from.size_; ++j) {
							to.indexToCoordinates(j, coords.begin());
							to(j) = from(coords.begin());
						}
					}
					else {
						for(size_t j=0; j<from.size_; ++j) {
							to.indexToCoordinates(j, coords.begin());
							to(coords.begin()) = from(coords.begin());
						}
					}
				}
			}
		}
	}
};

/// AssignmentOperatorHelper class
///
/// This is a helper class for the class template View. It is used
/// to implement a different behavoir of the assignment operator of
/// View, depending on the template parameter isConst of View.
/// Cf. the definition of operator= in View for details.
///
template<class T> 
struct AssignmentOperatorHelper<T, true>
{
	/// from either constant or mutable to constant
	template<bool isConstLocal>
	static void execute
	(
		const View<T, isConstLocal>& from,
		View<T, true>& to
	)
	{
		to.dimension_ = from.dimension_;
		to.shape_ = from.shape_;
		to.shapeStrides_ = from.shapeStrides_;
		to.size_ = from.size_;
		to.data_ = from.data_; // copy pointer
		to.offset_ = from.offset_;
		to.strides_ = from.strides_;
		to.isSimple_ = from.isSimple_; 
		to.coordinateOrder_ = from.coordinateOrder_;
	}
};
// \endcond end suppress doxygen

/// Assignment.
/// 
/// operator= (the assignment operator) has a non-trivial behavior.
/// In most cases, it will work as most programmers will expect.
/// Here's a complete description of the semantics of to.operator=(from)
/// or equivalently, to = from.
/// 
/// Consider the following cases:
/// (1) 'to' is mutable (isConst == false)
/// 	(a) 'from' is mutable (isConst == false)
/// 		(i) 'to' is initialized (data_ != 0)
/// 		(ii) 'to' is un-initialized (data_ == 0)
/// 	(b) 'from' is constant (isConst == true)
/// (2) 'to' is constant (isConst == true)
/// 
/// (i) The operator attempts to copy the data under view 'b' to
/// the memory under view 'a'. This works if both views have the
/// same size, regardless of their dimension and shape. Equality
/// of sizes is checked by an assertion.
/// 
/// (ii) Unless &a == &b (self-assignment), the operator copies
/// the (data) pointer of view 'b' to view 'a', without copying
/// the data itself. In addition, all the properties of view 'b'
/// are copied to view 'a'.
/// 
/// (b) The operator attempts to copy the data under view 'b' to
/// the memory under view 'a'. This works if both views have the
/// same size, regardless of their dimension and shape. Equality
/// of sizes is checked by an assertion. If 'a' is un-initialized
/// the assertion fails (because the size of a will be zero).
/// Unlike in (ii), the pointer is not copied in this case.
/// Thus, a conversion from mutable to const is prevented.
/// 
/// (2) Unless &a == &b (self-assignment), the operator copies
/// the (data) pointer of view 'b' to view 'a', without copying
/// the data itself. In addition, all the properties of view 'b'
/// are copied to view 'a'. Note that changing the data under
/// a constant view would be counter-intuitive.
/// 
template<class T, bool isConst> 
inline View<T, isConst>&
View<T, isConst>::operator=
(
	const View<T, true>& in
)
{
	testInvariant();
	if(this->overlaps(in)) {
		Marray<T> m = in; // temporary copy
		(*this) = m;
	}
	else {
		AssignmentOperatorHelper<T, isConst>::execute(in, *this);
	}
	testInvariant();
	return *this;
}

/// Assignment.
///
template<class T, bool isConst> 
inline View<T, isConst>&
View<T, isConst>::operator=
(
	const View<T, false>& in
)
{
	testInvariant();
	if(this->overlaps(in)) { 
		Marray<T> m = in; // temporary copy
		(*this) = m;
	}
	else {
		AssignmentOperatorHelper<T, isConst>::execute(in, *this);
	}
	testInvariant();
	return *this;
}

template<class T, bool isConst> 
inline View<T, isConst>& 
View<T, isConst>::operator=
(
	const T& value
)
{
	Assert(NO_DEBUG || data_ != 0);
	for(size_t j=0; j<size_; ++j) {
		(*this)(j) = value;
	}
	return *this;
}

/// Get a sub-view with the same coordinate order.
///
/// \param bit Iterator to the beginning of a coordinate sequence
/// that determines the start position of the sub-view.
/// \param sit Iterator to the beginning of a sequence
/// that determines the shape of the sub-view.
/// \param out Sub-View (output).
///
template<class T, bool isConst> 
template<class BaseIterator, class ShapeIterator>
inline void
View<T, isConst>::view
(
	BaseIterator bit,
	ShapeIterator sit,
	View<T, isConst>& out
) const
{
	view(bit, sit, coordinateOrder_, out);
}

/// Get a sub-view.
///
/// \param bit Iterator to the beginning of a coordinate sequence
/// that determines the start position of the sub-view.
/// \param sit Iterator to the beginning of a sequence
/// that determines the shape of the sub-view.
/// \param internalCoordinateOrder Flag to set the coordinate order
/// for scalar indexing and iterators of the sub-view.
/// \param out Sub-View (output).
///
template<class T, bool isConst> 
template<class BaseIterator, class ShapeIterator>
inline void
View<T, isConst>::view
(
	BaseIterator bit,
	ShapeIterator sit,
	const CoordinateOrder& internalCoordinateOrder,
	View<T, isConst>& out
) const
{
	testInvariant();
	size_t newOffset = 0;
	coordinatesToOffset(bit, newOffset);
	out.assign(sit, sit+dimension_, strides_.begin(),
		data_, newOffset, internalCoordinateOrder);
}

/// Get a sub-view with the same coordinate order.
///
/// \param bit Iterator to the beginning of a coordinate sequence
/// that determines the start position of the sub-view.
/// \param sit Iterator to the beginning of a sequence
/// that determines the shape of the sub-view.
/// \return Sub-View.
///
template<class T, bool isConst> 
template<class BaseIterator, class ShapeIterator>
inline View<T, isConst>
View<T, isConst>::view
(
    BaseIterator bit,
    ShapeIterator sit
) const
{
    View<T, isConst> v;
    this->view(bit, sit, v);
    return v;
}

/// Get a sub-view.
///
/// \param bit Iterator to the beginning of a coordinate sequence
/// that determines the start position of the sub-view.
/// \param sit Iterator to the beginning of a sequence
/// that determines the shape of the sub-view.
/// \param internalCoordinateOrder Flag to set the coordinate order
/// for scalar indexing and iterators of the sub-view.
/// \return Sub-View.
///
template<class T, bool isConst> 
template<class BaseIterator, class ShapeIterator>
inline View<T, isConst>
View<T, isConst>::view
(
    BaseIterator bit,
    ShapeIterator sit,
	const CoordinateOrder& internalCoordinateOrder
) const
{
    View<T, isConst> v;
    this->view(bit, sit, internalCoordinateOrder, v);
    return v;
}

/// Get a sub-view to constant data with the same coordinate
/// order.
///
/// \param bit Iterator to the beginning of a coordinate sequence
/// that determines the start position of the sub-view.
/// \param sit Iterator to the beginning of a sequence
/// that determines the shape of the sub-view.
/// \param out Sub-View (output).
///
template<class T, bool isConst> 
template<class BaseIterator, class ShapeIterator>
inline void
View<T, isConst>::constView
(
	BaseIterator bit,
	ShapeIterator sit,
	View<T, true>& out
) const
{
	constView(bit, sit, coordinateOrder_, out);
}

/// Get a sub-view to constant data.
///
/// \param bit Iterator to the beginning of a coordinate sequence
/// that determines the start position of the sub-view.
/// \param sit Iterator to the beginning of a sequence
/// that determines the shape of the sub-view.
/// \param internalCoordinateOrder Flag to set the coordinate order
/// for scalar indexing and iterators of the sub-view. 
/// \param out Sub-View (output).
///
template<class T, bool isConst> 
template<class BaseIterator, class ShapeIterator>
inline void
View<T, isConst>::constView
(
	BaseIterator bit,
	ShapeIterator sit,
	const CoordinateOrder& internalCoordinateOrder,
	View<T, true>& out
) const
{
	testInvariant();
	size_t newOffset = 0;
	coordinatesToOffset(bit, newOffset);
	out.assign(sit, sit+dimension_, strides_.begin(), (const T*)(data_),
		newOffset, internalCoordinateOrder);
}

/// Get a sub-view to constant data with the same coordinate
/// order.
///
/// \param bit Iterator to the beginning of a coordinate sequence
/// that determines the start position of the sub-view.
/// \param sit Iterator to the beginning of a sequence
/// that determines the shape of the sub-view.
/// \return Sub-View.
///
template<class T, bool isConst> 
template<class BaseIterator, class ShapeIterator>
inline View<T, true>
View<T, isConst>::constView
(
    BaseIterator bit,
    ShapeIterator sit
) const
{
    View<T, true> v;
    this->constView(bit, sit, v);
    return v;
}

/// Get a sub-view to constant data.
///
/// \param bit Iterator to the beginning of a coordinate sequence
/// that determines the start position of the sub-view.
/// \param sit Iterator to the beginning of a sequence
/// that determines the shape of the sub-view.
/// \param internalCoordinateOrder Flag to set the coordinate order
/// for scalar indexing and iterators of the sub-view. 
/// \return Sub-View.
///
template<class T, bool isConst> 
template<class BaseIterator, class ShapeIterator>
inline View<T, true>
View<T, isConst>::constView
(
    BaseIterator bit,
    ShapeIterator sit,
	const CoordinateOrder& internalCoordinateOrder
) const
{
    View<T, true> v;
    this->constView(bit, sit, internalCoordinateOrder, v);
    return v;
}

#ifdef HAVE_CPP0X_INITIALIZER_LISTS
/// Get a sub-view.
///
/// \param b Initializer list defining the coordinate sequence
/// that determines the start position of the sub-view.
/// \param s Initializer list defining the coordinate sequence
/// that determines the stop position of the sub-view.
/// \param internalCoordinateOrder Flag to set the coordinate order
/// for scalar indexing and iterators of the sub-view.
///
template<class T, bool isConst> 
inline void
View<T, isConst>::view
(
	std::initializer_list<size_t> b,
	std::initializer_list<size_t> s,
	const CoordinateOrder& internalCoordinateOrder,
	View<T, isConst>& out
) const
{
	view(b.begin(), s.begin(), internalCoordinateOrder, out);
}

/// Get a sub-view with the same coordinate order.
///
/// \param b Initializer list coordinate sequence
/// that determines the start position of the sub-view.
/// \param s Initializer list coordinate sequence
/// that determines the stop position of the sub-view.
/// \param out Sub-View (output).
///
template<class T, bool isConst> 
inline void
View<T, isConst>::view
(
	std::initializer_list<size_t> b,
	std::initializer_list<size_t> s,
	View<T, isConst>& out
) const 
{
	view(b.begin(), s.begin(), coordinateOrder_, out);
}

/// Get a sub-view to constant data.
///
/// \param b Initializer list coordinate sequence
/// that determines the start position of the sub-view.
/// \param s Initializer list coordinate sequence
/// that determines the stop position of the sub-view.
/// \param internalCoordinateOrder Flag to set the coordinate order
/// for scalar indexing and iterators of the sub-view. 
///
template<class T, bool isConst> 
inline void
View<T, isConst>::constView
(
	std::initializer_list<size_t> b,
	std::initializer_list<size_t> s,
	const CoordinateOrder& internalCoordinateOrder,
	View<T, true>& out
) const
{
	constView(b.begin(), s.begin(), internalCoordinateOrder, out);
}

/// Get a sub-view to constant data with the same coordinate
/// order.
///
/// \param b Initializer list coordinate sequence
/// that determines the start position of the sub-view.
/// \param s Initializer list coordinate sequence
/// that determines the stop position of the sub-view.
/// \param out Sub-View (output).
///
template<class T, bool isConst>
inline void
View<T, isConst>::constView
(
	std::initializer_list<size_t> b,
	std::initializer_list<size_t> s,
	View<T, true>& out
) const
{
	constView(b.begin(), s.begin(), coordinateOrder_, out);
}
#endif

/// Reshape the View.
/// 
/// Two conditions have to be fulfilled in order for reshape to work:
/// - The new and the old shape must have the same size.
/// - The view must be simple, cf. isSimple().
/// .
/// 
/// \param begin Iterator to the beginning of a sequence that determines
/// the new shape.
/// \param end Iterator to the end of that sequence.
///
/// \sa reshapedView(), isSimple()
///
template<class T, bool isConst> 
template<class ShapeIterator>
inline void
View<T, isConst>::reshape
(
	ShapeIterator begin,
	ShapeIterator end
)
{
	testInvariant();
	Assert(NO_DEBUG || isSimple());
		// only simple views can be reshaped! this includes all
		// Marrays, Matrices, and Vectors.
	if(!NO_ARG_TEST) {
		ShapeIterator it = begin;
		size_t size = 1;
		while(it != end) {
			Assert(static_cast<size_t>(*it) > 0);
			size *= static_cast<size_t>(*it);
			++it;
		}
		Assert(size == size_);
	}
	assign(begin, end, data_, coordinateOrder_, coordinateOrder_);
	testInvariant();
}

/// Get a reshaped View.
/// 
/// Two conditions have to be fulfilled:
/// - The new and the old shape must have the same size.
/// - The view must be simple, cf. isSimple().
/// .
/// 
/// \param begin Iterator to the beginning of a sequence that determines
/// the new shape.
/// \param end Iterator to the end of that sequence.
///
/// \sa reshape(), isSimple()
///
template<class T, bool isConst> 
template<class ShapeIterator>
inline View<T, isConst>
View<T, isConst>::reshapedView
(
	ShapeIterator begin,
	ShapeIterator end
) const
{
	View<T, isConst> out = *this;
	out.reshape(begin, end);
	return out;
}

#ifdef HAVE_CPP0X_INITIALIZER_LISTS
/// Reshape the View.
/// 
/// Two conditions have to be fulfilled in order for reshape to work:
/// - The new and the old shape must have the same size.
/// - The view must be simple, cf. isSimple().
/// .
/// 
/// \param shape Initializer list defining the new shape.
///
/// \sa reshapedView(), isSimple()
///
template<class T, bool isConst> 
inline void
View<T, isConst>::reshape
(
	std::initializer_list<size_t> shape
)
{
	reshape(shape.begin(), shape.end());
}

/// Get a reshaped View.
/// 
/// Two conditions have to be fulfilled:
/// - The new and the old shape must have the same size.
/// - The view must be simple, cf. isSimple().
/// .
/// 
/// \param shape Initializer list defining the new shape.
///
/// \sa reshape(), isSimple()
///
template<class T, bool isConst> 
inline View<T, isConst>
View<T, isConst>::reshapedView
(
	std::initializer_list<size_t> shape
) const
{
	return reshapedView(shape.begin(), shape.end());
}
#endif

/// Get a View where one coordinate is bound to a value.
///
/// Binds one coordinate to a certain value. This reduces the
/// dimension by 1.
///
/// \param dimension Dimension of the coordinate to bind.
/// \param value Value to assign to the coordinate.
/// \return The bound view.
/// \sa squeeze(), squeezeView()
///
template<class T, bool isConst> 
View<T, isConst>
View<T, isConst>::boundView
(
	const dimension_type& dimension,
	const coordinate_type& value
) const
{
	testInvariant();
	Assert(NO_ARG_TEST || (dimension < dimension_
		&& value < shape_[static_cast<size_t>(dimension)]));

	View v = *this; // copy
	{
		coordinate_tuple coord(dimension_);
		coord[dimension] = value;
		size_t newOffset;
		coordinatesToOffset(coord.begin(), newOffset);
		v.offset_ = newOffset;
	}
	v.size_ /= v.shape_[dimension];
	v.shape_.erase(v.shape_.begin() + dimension);
	if(v.dimension_ == 1) {
		v.shapeStrides_.erase(v.shapeStrides_.begin());
	}
	else {
		v.shapeStrides_.erase(v.shapeStrides_.end() - 1);
		shapeStrides(v.shape_.begin(), v.shape_.end(),
			v.shapeStrides_.begin(), v.coordinateOrder_);
	}
	v.strides_.erase(v.strides_.begin() + dimension);
	--v.dimension_;
	v.updateSimplicity();
	v.testInvariant();
	return v;
}

/// Remove singleton dimensions by setting their coordinates to zero.
///
/// \sa squeezedView(), boundView()
///
template<class T, bool isConst> 
void
View<T, isConst>::squeeze()
{
	testInvariant();
	dimension_type j = 0;
	while(j < dimension_) {
		if(shape_[static_cast<size_t>(j)] == 1) {
			shape_.erase(shape_.begin() + j);
			strides_.erase(strides_.begin() + j);
			--dimension_;
		}
		else {
			++j;
		}
	}
	shapeStrides_.resize(dimension_);
	if(dimension_ != 0) {
		shapeStrides(shape_.begin(), shape_.end(), shapeStrides_.begin(),
			coordinateOrder_);
	}
	updateSimplicity();
	testInvariant();
}

/// Get a View where all singleton dimensions have been removed by
/// setting their coordinates to zero.
///
/// \sa squeeze(), boundView()
///
template<class T, bool isConst> 
inline View<T, isConst>
View<T, isConst>::squeezedView() const
{
	View<T, isConst> v = *this;
	v.squeeze();
	return v;
}

#ifdef HAVE_CPP0X_INITIALIZER_LISTS
/// Permute dimensions.
///
/// \param begin Iterator to the beginning of a sequence which
/// has to contain the integers 0, ..., dimension()-1 in any
/// order. Otherwise, a runtime error is thrown.
/// \sa permutedView(), transpose(), transposedView(), shift(),
/// shiftedView()
///
template<class T, bool isConst>
void
View<T, isConst>::permute
(
	std::initializer_list<size_t> permutation
)
{
	permute(permutation.begin());
}
#endif

/// Permute dimensions.
///
/// \param begin Iterator to the beginning of a sequence which
/// has to contain the integers 0, ..., dimension()-1 in any
/// order. Otherwise, a runtime error is thrown.
/// \sa permutedView(), transpose(), transposedView(), shift(),
/// shiftedView()
///
template<class T, bool isConst> 
template<class CoordinateIterator>
void
View<T, isConst>::permute
(
	CoordinateIterator begin
)
{
	testInvariant();
	if(!NO_ARG_TEST) {
		Assert(dimension_ != 0);
		std::set<coordinate_type> s1, s2;
		CoordinateIterator it = begin;
		for(size_t j=0; j<dimension_; ++j) {
			s1.insert(j);
			s2.insert(*it);
			++it;
		}
		Assert(s1 == s2);
	}

	// compute new shape and new strides
	coordinate_tuple newShape = coordinate_tuple(dimension_);
	offset_tuple newStrides = offset_tuple(dimension_);
	for(size_t j=0; j<dimension_; ++j) {
		newShape[j] = shape_[static_cast<size_t>(*begin)];
		newStrides[j] = strides_[static_cast<size_t>(*begin)];
		++begin;
	}

	// update shape, shape strides, strides, and simplicity
	shape_ = newShape;
	shapeStrides(shape_.begin(), shape_.end(),
		shapeStrides_.begin(), coordinateOrder_);
	strides_ = newStrides;
	updateSimplicity();

	testInvariant();
}

/// Get a View with permuted dimensions.
///
/// \param begin Iterator to the beginning of a sequence which
/// has to contain the integers 0, ..., dimension()-1 in any
/// order. Otherwise, a runtime error is thrown.
/// \return Permuted View.
/// \sa permute(), transpose(), transposedView(), shift(),
/// shiftedView()
///
template<class T, bool isConst> 
template<class CoordinateIterator>
inline View<T, isConst>
View<T, isConst>::permutedView
(
	CoordinateIterator begin
) const
{
	View<T, isConst> out = *this;
	out.permute(begin);
	return out;
}

/// Exchange two dimensions.
///
/// \param c1 Dimension
/// \param c2 Dimension
/// \sa permute(), permutedView(), shift(), shiftedView()
///
template<class T, bool isConst> 
void
View<T, isConst>::transpose
(
	const coordinate_type& c1,
	const coordinate_type& c2
)
{
	testInvariant();
	Assert(NO_ARG_TEST ||
		(dimension_ != 0 && c1 < dimension_ && c2 < dimension_));

	size_t j1 = static_cast<size_t>(c1);
	size_t j2 = static_cast<size_t>(c2);
	coordinate_type c;
	size_t d;

	// transpose shape
	c = shape_[j2];
	shape_[j2] = shape_[j1];
	shape_[j1] = c;

	// transpose strides
	d = strides_[j2];
	strides_[j2] = strides_[j1];
	strides_[j1] = d;

	// update shape strides
	shapeStrides(shape_.begin(), shape_.end(),
		shapeStrides_.begin(), coordinateOrder_);

	updateSimplicity();
	testInvariant();
}

/// Reverse dimensions.
///
/// \sa transposedView(), permute(), permutedView(), shift(),
/// shiftedView()
///
template<class T, bool isConst> 
void
View<T, isConst>::transpose()
{
	testInvariant();
	for(size_t j=0; j<static_cast<size_t>(dimension_)/2; ++j) {
		size_t k = static_cast<size_t>(dimension_)-j-1;

		// transpose shape
		size_t tmp = shape_[j];
		shape_[j] = shape_[k];
		shape_[k] = tmp;

		// transpose strides
		tmp = strides_[j];
		strides_[j] = strides_[k];
		strides_[k] = tmp;
	}
	shapeStrides(shape_.begin(), shape_.end(),
		shapeStrides_.begin(), coordinateOrder_);
	updateSimplicity();
	testInvariant();
}

/// Get a View with two dimensions exchanged.
///
/// \param c1 Dimension
/// \param c2 Dimension
/// \return Transposed View.
/// \sa transpose(), permute(), permutedView(), shift(),
/// shiftedView()
///
template<class T, bool isConst> 
inline View<T, isConst>
View<T, isConst>::transposedView
(
	const coordinate_type& c1,
	const coordinate_type& c2
) const
{
	View<T, isConst> out = *this;
	out.transpose(c1, c2);
	return out;
}

/// Get a View with dimensions reversed.
///
/// \return View with dimensions reversed.
/// \sa transpose(), permute(), permutedView(), shift(),
/// shiftedView()
///
template<class T, bool isConst> 
inline View<T, isConst>
View<T, isConst>::transposedView() const
{
	View<T, isConst> out = *this;
	out.transpose();
	return out;
}

/// Cycle shift dimensions.
///
/// \param n Number of positions to shift
/// \sa shiftedView(), permute(), permutedView(), transpose(),
/// transposedView()
///
template<class T, bool isConst> 
inline void
View<T, isConst>::shift
(
	int n
) 
{
	testInvariant();
	Assert(NO_DEBUG || dimension_ != 0);

	if(n <= -static_cast<int>(dimension_) || n >= static_cast<int>(dimension_)) {
		shift(n % static_cast<int>(dimension_));
	}
	else {
		if(n > 0) {
			shift(n - static_cast<int>(dimension_));
		}
		else {
			coordinate_tuple p(dimension_);
			for(size_t j=0; j<dimension_; ++j) {
				p[j] = static_cast<size_t>((static_cast<int>(j) - n)) % dimension_;
			}
			permute(p.begin());
		}
	}

	testInvariant();
}

/// Get a View which dimensions cycle shifted.
///
/// \param n Number of positions to shift
/// \sa shift(), permute(), permutedView(), transpose(), transposedView()
///
template<class T, bool isConst> 
inline View<T, isConst>
View<T, isConst>::shiftedView
(
	int n
) const
{
	View<T, isConst> out = *this;
	out.shift(n);
	return out;
}

#ifdef MARRAY_COMPATIBILITY
/// For compatiblity with Vigra MultiArray.
///
/// \param bit Iterator to the beginning of a coordinate sequence
/// that determines the start position of the sub-view.
/// \param eit Iterator to the beginning of a coordinate sequence
/// that determines the stop position of the sub-view PLUS 1 in
/// each dimension.
/// \return View ranging from bit to eit-1
///
template<class T, bool isConst> 
template<class BaseIterator, class EndIterator>
View<T, isConst>
View<T, isConst>::subarray
(
	BaseIterator bit,
	EndIterator eit
) const
{
	coordinate_tuple base(dimension_);
	coordinate_tuple shape(dimension_);
	for(size_t j=0; j<dimension_; ++j, ++bit, ++eit) {
		Assert(NO_ARG_TEST || *bit < *eit);
		base[j] = *bit;
		shape[j] = *eit - *bit;
	}
	return this->view(base.begin(), shape.begin());
}
#endif // #ifdef MARRAY_COMPATIBILITY

/// Get an iterator to the beginning.
///
/// \return Iterator.
/// \sa end()
///
template<class T, bool isConst> 
inline typename View<T, isConst>::iterator
View<T, isConst>::begin()
{
	testInvariant();
	return Iterator<T, isConst>(*this, 0);
}

/// Get the end-iterator.
///
/// \return Iterator.
/// \sa begin()
///
template<class T, bool isConst> 
inline typename View<T, isConst>::iterator
View<T, isConst>::end()
{
	testInvariant();
	return Iterator<T, isConst>(*this, size_);
}

/// Get an iterator to the beginning.
///
/// \return Iterator.
/// \sa end()
///
template<class T, bool isConst> 
inline typename View<T, isConst>::const_iterator
View<T, isConst>::begin() const
{
	testInvariant();
	return Iterator<T, true>(*this, 0);
}

/// Get the end-iterator.
///
/// \return Iterator.
/// \sa begin()
///
template<class T, bool isConst> 
inline typename View<T, isConst>::const_iterator
View<T, isConst>::end() const
{
	testInvariant();
	return Iterator<T, true>(*this, size_);
}

/// Get a reserve iterator to the beginning.
///
/// \return Iterator.
/// \sa rend()
///
template<class T, bool isConst> 
inline typename View<T, isConst>::reverse_iterator
View<T, isConst>::rbegin()
{
	return reverse_iterator(end());
}

/// Get the reverse end-iterator.
///
/// \return Iterator.
/// \sa rbegin()
///
template<class T, bool isConst> 
inline typename View<T, isConst>::reverse_iterator
View<T, isConst>::rend()
{
	return reverse_iterator(begin());
}

/// Get a reserve iterator to the beginning.
///
/// \return Iterator.
/// \sa rend()
///
template<class T, bool isConst> 
inline typename View<T, isConst>::const_reverse_iterator
View<T, isConst>::rbegin() const
{
	return const_reverse_iterator(end());
}

/// Get the reverse end-iterator.
///
/// \return Iterator.
/// \sa rbegin()
///
template<class T, bool isConst> 
inline typename View<T, isConst>::const_reverse_iterator
View<T, isConst>::rend() const
{
	return const_reverse_iterator(begin());
}

/// Update Simplicity.
///
/// This function sets the redundant boolean attribute isSimple_.
/// isSimple_ is set to true if the view is simple, i.e. if the
/// shape strides equal the strides and the offset is zero. 
///
template<class T, bool isConst> 
void
View<T, isConst>::updateSimplicity()
{
	// the invariant is not tested in this function
	// because this function is called during the
	// construction of a view
	if(offset_ == 0) {
		isSimple_ = true;
		for(size_t j=0; j<dimension_; ++j) {
			if(shapeStrides_[j] != strides_[j]) {
				isSimple_ = false;
				break;
			}
		}
	}
	else {
		isSimple_ = false;
	}
}

/// Test invariant.
///
/// This function tests the invariant of View and thus the consistency
/// of redundant information.
///
template<class T, bool isConst> 
inline void
View<T, isConst>::testInvariant() const
{
	if(!NO_DEBUG) {
		if(dimension_ == 0) {
			if(data_ == 0) {
				// un-initialized
				Assert(shape_.size() == 0 && shapeStrides_.size() == 0
					&& size_ == 0 && offset_ == 0 && strides_.size() == 0 
					&& isSimple_ == true);
			}
			else {
				// scalar
				Assert(shape_.size() == 0 && shapeStrides_.size() == 0
					&& size_ == 1 && strides_.size() == 0);
			}
		}
		else {
			Assert(data_ != 0 && shape_.size() == static_cast<size_t>(dimension_)
				&& shapeStrides_.size() == static_cast<size_t>(dimension_)
				&& strides_.size() == static_cast<size_t>(dimension_));

			// test size_ to be consistent with shape_
			size_t testSize = 1;
			for(size_t j=0; j<static_cast<size_t>(dimension_); ++j) {
				testSize *= shape_[j];
			}
			Assert(size_ == testSize);

			// test shapeStrides_ to be consistent with shape_
			if(coordinateOrder_ == FirstMajorOrder) {
				size_t tmp = 1;
				for(size_t j=0; j<dimension_; ++j) {
					Assert(shapeStrides_[dimension_-j-1] == tmp);
					tmp *= shape_[dimension_-j-1];
				}
			}
			else {
				size_t tmp = 1;
				for(size_t j=0; j<dimension_; ++j) {
					Assert(shapeStrides_[j] == tmp);
					tmp *= shape_[j];
				}
			}

			// test the simplicity condition 
			if(isSimple_) {
				Assert(offset_ == 0);
				for(size_t j=0; j<static_cast<size_t>(dimension_); ++j) {
					Assert(strides_[j] == shapeStrides_[j]);
				}
			}
		}
	}
}

/// Check whether two Views overlap.
///
/// This function returns true if two memory intervals overlap:
/// (1) the interval between the first and the last element of the object
/// whose member function overlaps() is called.
/// (2) the interval between the first and the last element of v.
///
/// Note that this not necessarily implies the existence of an element 
/// that is addressed by both v and the current object. v could for
/// instance address all odd elements in a vector while the current
/// object addresses all even elements. 
///
/// \param v A view to compare with *this.
/// \return bool.
///
template<class T, bool isConst> 
template<bool isConstLocal>
inline bool View<T, isConst>::overlaps
(
	const View<T, isConstLocal>& v
)
{
	testInvariant();
	if(!NO_ARG_TEST) {
		v.testInvariant();
	}

	if(data_ == 0 || v.data_ == 0) {
		return false;
	}
	else {
		const T* maxPointer = & (*this)(size_-1);
		const T* maxPointerV = & v(v.size_-1);
		if(    (data_   <= v.data_ && v.data_ <= maxPointer)
		    || (v.data_ <= data_   && data_   <= maxPointerV) )
		{
			return true;
		}
	}
	return false;
}

/// Output as string.
///
template<class T, bool isConst>
std::string 
View<T, isConst>::asString
(
	const StringStyle& style
) const
{
	testInvariant();
	std::ostringstream out(std::ostringstream::out);
	if(style == MatrixStyle) {
		if(dimension_ == 0) {
			// scalar
			out << "A = " << (*this)(0) << std::endl;
		}
		else if(dimension_ == 1) {
			// vector
			out << "A = (";
			for(size_t j=0; j<this->size(); ++j) {
				out << (*this)(j) << ", ";
			}
			out << "\b\b)" << std::endl;
		}
		else if(dimension_ == 2) {
			// matrix
			if(coordinateOrder_ == FirstMajorOrder) {
				out << "A(r,c) =" << std::endl;
				for(size_t y=0; y<this->shape(0); ++y) {
					for(size_t x=0; x<this->shape(1); ++x) {
						out << (*this)(y, x) << ' ';
					}
					out << std::endl;
				}
			}
			else {
				out << "A(c,r) =" << std::endl;
				for(size_t y=0; y<this->shape(1); ++y) {
					for(size_t x=0; x<this->shape(0); ++x) {
						out << (*this)(x, y) << ' ';
					}
					out << std::endl;
				}
			}
		}
		else {
			// higher dimensional
			std::vector<coordinate_type> c1(dimension_);
			std::vector<coordinate_type> c2(dimension_);
			unsigned short q = 2;
			if(coordinateOrder_ == FirstMajorOrder) {
				q = dimension_-3;
			}
			for(const_iterator it = this->begin(); it.hasMore(); ++it) {
				it.coordinate(c2.begin());
				if(it.index() == 0 || c2[q] != c1[q]) {
					if(it.index() != 0) {
						out << std::endl << std::endl;
					}
					if(coordinateOrder_ == FirstMajorOrder) {
						out << "A(";
						for(size_t j=0; j<static_cast<size_t>(dimension_-2); ++j) {
							out << c2[j] << ",";
						}
					}
					else {
						out << "A(c,r,";
						for(size_t j=2; j<dimension_; ++j) {
							out << c2[j] << ",";
						}
					}
					out << '\b';
					if(coordinateOrder_ == FirstMajorOrder) {
						out << ",r,c";
					}
					out << ") =" << std::endl;
				}
				else if(c2[1] != c1[1]) {
					out << std::endl;
				}
				out << *it << " ";
				c1 = c2;
			}
			out << std::endl;
		}
		out << std::endl;
	}
	else if(style == TableStyle) {
		if(dimension_ == 0) {
			// scalar
			out << "A = " << (*this)(0) << std::endl;
		}
		else {
			// non-scalar
			std::vector<coordinate_type> c(dimension_);
			for(const_iterator it = this->begin(); it.hasMore(); ++it) {
				out << "A(";
				it.coordinate(c.begin());
				for(size_t j=0; j<c.size(); ++j) {
					out << c[j] << ',';
				}
				out << "\b) = " << *it << std::endl;
			}
		}
		out << std::endl;
	}
	return out.str();
}

#ifdef MARRAY_COMPATIBILITY
/// For compatiblity with Vigra MultiArray, equivalent to permute()
///
template<class T, bool isConst>
template<class CoordinateIterator>
inline View<T, isConst>
View<T, isConst>::permuteDimensions
(
	CoordinateIterator it
) const
{
	return this->permutedView(it);
}

/// For compatiblity with Vigra MultiArray, equivalent to shift() 
///
template<class T, bool isConst>
inline View<T, isConst>
View<T, isConst>::shiftDimensions
(
	int z
) const
{
	return this->shiftedView(z);
}
#endif // #ifdef MARRAY_COMPATIBILITY

// implementation of arithmetic operators of View

template<class T>
inline View<T, false>&
operator+=
(
	View<T, false>& v,
	const T& x
)
{
	for(size_t j=0; j<v.size(); ++j) {
		v(j) += x;
	}
	return v;
}

template<class T1, class T2, bool isConst>
inline View<T1, false>&
operator+=
(
	View<T1, false>& v,
	const View<T2, isConst>& w
)
{
	if(!NO_ARG_TEST) {
		Assert(w.dimension() == 0 || v.dimension() == w.dimension());
		if(w.dimension() != 0) {
			for(size_t j=0; j<v.dimension(); ++j) {
				Assert(v.shape(j) == w.shape(j));
			}
		}
	}
	if(w.dimension() == 0) {
		for(size_t j=0; j<v.size(); ++j) {
			v(j) += w(0);
		}
	}
	else {
		if(v.overlaps(w)) {
			Marray<T2> m = w; // temporary copy
			v += m; 
		}
		else {
			if(v.coordinateOrder() == w.coordinateOrder()) {
				for(size_t j=0; j<v.size(); ++j) {
					v(j) += w(j);
				}
			}
			else {
				std::vector<size_t> coords(v.dimension());
				if(v.isSimple()) {
					for(size_t j=0; j<v.size(); ++j) {
						v.indexToCoordinates(j, coords.begin());
						v(j) += w(coords.begin());
					}
				}
				else if(w.isSimple()) {
					for(size_t j=0; j<v.size(); ++j) {
						w.indexToCoordinates(j, coords.begin());
						v(coords.begin()) += w(j);
					}
				}
				else {
					for(size_t j=0; j<v.size(); ++j) {
						v.indexToCoordinates(j, coords.begin());
						v(coords.begin()) += w(coords.begin());
					}
				}
			}
		}
	}
	return v;
}

// prefix
template<class T>
inline View<T, false>&
operator++
(
	View<T, false>& v
)
{
	for(size_t j=0; j<v.size(); ++j) {
		++v(j);
	}
	return v;
}

// postfix
template<class T>
inline Marray<T>
operator++
(
	Marray<T>& in,
	int dummy
) 
{
	Marray<T> out = in; 
	for(size_t j=0; j<in.size(); ++j) {
		++in(j);
	}
	return out;
}

template<class T>
inline View<T, false>&
operator-=
(
	View<T, false>& v,
	const T& x
)
{
	for(size_t j=0; j<v.size(); ++j) {
		v(j) -= x;
	}
	return v;
}

template<class T1, class T2, bool isConst>
inline View<T1, false>&
operator-=
(
	View<T1, false>& v,
	const View<T2, isConst>& w
)
{
	if(!NO_ARG_TEST) {
		Assert(w.dimension() == 0 || v.dimension() == w.dimension());
		if(w.dimension() != 0) {
			for(size_t j=0; j<v.dimension(); ++j) {
				Assert(v.shape(j) == w.shape(j));
			}
		}
	}
	if(w.dimension() == 0) {
		for(size_t j=0; j<v.size(); ++j) {
			v(j) -= w(0);
		}
	}
	else {
		if(v.overlaps(w)) {
			Marray<T2> m = w; // temporary copy
			v -= m; 
		}
		else {
			if(v.coordinateOrder() == w.coordinateOrder()) {
				for(size_t j=0; j<v.size(); ++j) {
					v(j) -= w(j);
				}
			}
			else {
				std::vector<size_t> coords(v.dimension());
				if(v.isSimple()) {
					for(size_t j=0; j<v.size(); ++j) {
						v.indexToCoordinates(j, coords.begin());
						v(j) -= w(coords.begin());
					}
				}
				else if(w.isSimple()) {
					for(size_t j=0; j<v.size(); ++j) {
						w.indexToCoordinates(j, coords.begin());
						v(coords.begin()) -= w(j);
					}
				}
				else {
					for(size_t j=0; j<v.size(); ++j) {
						v.indexToCoordinates(j, coords.begin());
						v(coords.begin()) -= w(coords.begin());
					}
				}
			}
		}
	}
	return v;
}

// prefix
template<class T>
inline View<T, false>&
operator--
(
	View<T, false>& v
)
{
	for(size_t j=0; j<v.size(); ++j) {
		--v(j);
	}
	return v;
}

// postfix
template<class T>
inline Marray<T>
operator--
(
	Marray<T>& in,
	int dummy
) 
{
	Marray<T> out = in; 
	for(size_t j=0; j<in.size(); ++j) {
		--in(j);
	}
	return out;
}

template<class T>
inline View<T, false>&
operator*=
(
	View<T, false>& v,
	const T& x
)
{
	for(size_t j=0; j<v.size(); ++j) {
		v(j) *= x;
	}
	return v;
}

template<class T1, class T2, bool isConst>
inline View<T1, false>&
operator*=
(
	View<T1, false>& v,
	const View<T2, isConst>& w
)
{
	if(!NO_ARG_TEST) {
		Assert(w.dimension() == 0 || v.dimension() == w.dimension());
		if(w.dimension() != 0) {
			for(size_t j=0; j<v.dimension(); ++j) {
				Assert(v.shape(j) == w.shape(j));
			}
		}
	}
	if(w.dimension() == 0) {
		for(size_t j=0; j<v.size(); ++j) {
			v(j) *= w(0);
		}
	}
	else {
		if(v.overlaps(w)) {
			Marray<T2> m = w; // temporary copy
			v *= m; 
		}
		else {
			if(v.coordinateOrder() == w.coordinateOrder()) {
				for(size_t j=0; j<v.size(); ++j) {
					v(j) *= w(j);
				}
			}
			else {
				std::vector<size_t> coords(v.dimension());
				if(v.isSimple()) {
					for(size_t j=0; j<v.size(); ++j) {
						v.indexToCoordinates(j, coords.begin());
						v(j) *= w(coords.begin());
					}
				}
				else if(w.isSimple()) {
					for(size_t j=0; j<v.size(); ++j) {
						w.indexToCoordinates(j, coords.begin());
						v(coords.begin()) *= w(j);
					}
				}
				else {
					for(size_t j=0; j<v.size(); ++j) {
						v.indexToCoordinates(j, coords.begin());
						v(coords.begin()) *= w(coords.begin());
					}
				}
			}
		}
	}
	return v;
}

template<class T>
inline View<T, false>&
operator/=
(
	View<T, false>& v,
	const T& x
)
{
	for(size_t j=0; j<v.size(); ++j) {
		v(j) /= x;
	}
	return v;
}

template<class T1, class T2, bool isConst>
inline View<T1, false>&
operator/=
(
	View<T1, false>& v,
	const View<T2, isConst>& w
)
{
	if(!NO_ARG_TEST) {
		Assert(w.dimension() == 0 || v.dimension() == w.dimension());
		if(w.dimension() != 0) {
			for(size_t j=0; j<v.dimension(); ++j) {
				Assert(v.shape(j) == w.shape(j));
			}
		}
	}
	if(w.dimension() == 0) {
		for(size_t j=0; j<v.size(); ++j) {
			v(j) /= w(0);
		}
	}
	else {
		if(v.overlaps(w)) {
			Marray<T2> m = w; // temporary copy
			v /= m; 
		}
		else {
			if(v.coordinateOrder() == w.coordinateOrder()) {
				for(size_t j=0; j<v.size(); ++j) {
					v(j) /= w(j);
				}
			}
			else {
				std::vector<size_t> coords(v.dimension());
				if(v.isSimple()) {
					for(size_t j=0; j<v.size(); ++j) {
						v.indexToCoordinates(j, coords.begin());
						v(j) /= w(coords.begin());
					}
				}
				else if(w.isSimple()) {
					for(size_t j=0; j<v.size(); ++j) {
						w.indexToCoordinates(j, coords.begin());
						v(coords.begin()) /= w(j);
					}
				}
				else {
					for(size_t j=0; j<v.size(); ++j) {
						v.indexToCoordinates(j, coords.begin());
						v(coords.begin()) /= w(coords.begin());
					}
				}
			}
		}
	}
	return v;
}

template<class T, bool isConst>
inline Marray<T>
operator+
(
	const View<T, isConst>& v,
	const T& x
)
{
	Marray<T> m = v;
	m += x;
	return m;
}

template<class T, bool isConst>
inline Marray<T>
operator+
(
	const T& x,
	const View<T, isConst>& v
)
{
	Marray<T> m = v;
	m += x;
	return m;
}

template<class T, bool isConst1, bool isConst2>
inline Marray<T>
operator+
(
	const View<T, isConst1>& v,
	const View<T, isConst2>& w
)
{
	if(v.dimension() == 0) {
		Marray<T> m = w;
		m += v;
		return m;
	}
	else {
		Marray<T> m = v;
		m += w;
		return m;
	}
}

// unary
template<class T, bool isConst>
inline Marray<T>
operator+
(
	const View<T, isConst>& v
) 
{
	Marray<T> m = v;
	return m;
}

template<class T, bool isConst>
inline Marray<T>
operator-
(
	const View<T, isConst>& v,
	const T& x
)
{
	Marray<T> m = v;
	m -= x;
	return m;
}

template<class T, bool isConst>
inline Marray<T>
operator-
(
	const T& x,
	const View<T, isConst>& v
)
{
	Marray<T> m = v;
    for(size_t j=0; j<v.size(); ++j) {
        m(j) = x - v(j);
    }
	return m;
}

template<class T, bool isConst1, bool isConst2>
inline Marray<T>
operator-
(
	const View<T, isConst1>& v,
	const View<T, isConst2>& w
)
{
	if(v.dimension() == 0) {
        if(w.dimension() == 0) {
		    Marray<T> m = v;
		    m -= w(0);
		    return m;
        }
        else {
		    Marray<T> m = w;
            for(size_t j=0; j<w.size(); ++j) {
                m(j) = v(0) - w(j); 
            }
            return m;
        }
	}
	else {
		Marray<T> m = v;
		m -= w;
		return m;
	}
}

// unary
template<class T, bool isConst>
inline Marray<T>
operator-
(
	const View<T, isConst>& v
) 
{
	Marray<T> m = v;
	for(size_t j=0; j<m.size(); ++j) {
		m(j) = -m(j);
	}
	return m;
}

template<class T, bool isConst>
inline Marray<T>
operator*
(
	const View<T, isConst>& v,
	const T& x
)
{
	Marray<T> m = v;
	m *= x;
	return m;
}

template<class T, bool isConst>
inline Marray<T>
operator*
(
	const T& x,
	const View<T, isConst>& v
)
{
	Marray<T> m = v;
	m *= x;
	return m;
}

template<class T, bool isConst1, bool isConst2>
inline Marray<T>
operator*
(
	const View<T, isConst1>& v,
	const View<T, isConst2>& w
)
{
	if(v.dimension() == 0) {
		Marray<T> m = w;
		m *= v;
		return m;
	}
	else {
		Marray<T> m = v;
		m *= w;
		return m;
	}
}

template<class T, bool isConst>
inline Marray<T>
operator/
(
	const View<T, isConst>& v,
	const T& x
)
{
	Marray<T> m = v;
	m /= x;
	return m;
}

template<class T, bool isConst>
inline Marray<T>
operator/
(
	const T& x,
	const View<T, isConst>& v
)
{
	Marray<T> m = v;
    for(size_t j=0; j<v.size(); ++j) {
        m(j) = x / v(j);
    }
	return m;
}

template<class T, bool isConst1, bool isConst2>
inline Marray<T>
operator/
(
	const View<T, isConst1>& v,
	const View<T, isConst2>& w
)
{
	if(v.dimension() == 0) {
        if(w.dimension() == 0) {
		    Marray<T> m = v;
		    m /= w(0);
		    return m;
        }
        else {
		    Marray<T> m = w;
            for(size_t j=0; j<w.size(); ++j) {
                m(j) = v(0) / w(j); 
            }
            return m;
        }
	}
	else {
		Marray<T> m = v;
		m /= w;
		return m;
	}
}

// implementation of Marray

/// Clear Marray.
///
/// Leaves the Marray in the same state as if the empty constructor
/// had been called. Previously allocated memory is de-allocated.
///
/// \sa Marray()
///
template<class T> 
inline void
Marray<T>::assign()
{
  if(this->data_ != 0) {
    delete[] this->data_;
    this->data_ = 0;
  }
  base::assign();
}

/// Empty constructor.
///
template<class T> 
inline
Marray<T>::Marray()
{
	testInvariant();
}

/// Construct 0-dimensional (scalar) array.
///
/// \param value Value of the single data item.
/// \param coordinateOrder Flag specifying whether FirstMajorOrder or
/// LastMajorOrder is to be used. As the Marray can be resized after 
/// construction, the coordinate order has to be set even for a
/// 0-dimensional Marray.
///
template<class T> 
inline
Marray<T>::Marray
(
	const T& value,
	const CoordinateOrder& coordinateOrder
)
{
	this->data_ = new value_type[1];
	this->data_[0] = value;
	this->size_ = 1;
	this->coordinateOrder_ = coordinateOrder;
	testInvariant();
}

/// Copy from a Marray.
///
/// \param in Marray (source).
///
template<class T> 
Marray<T>::Marray
(
	const Marray<T>& in
)
{
	if(!NO_ARG_TEST) {
		in.testInvariant();
	}

	this->dimension_ = in.dimension_;
	this->shape_ = in.shape_;
	this->shapeStrides_ = in.shapeStrides_;
	this->size_ = in.size_;
	this->offset_ = 0;
	this->strides_ = in.strides_;
	this->isSimple_ = true;
	this->coordinateOrder_ = in.coordinateOrder_;
	if(in.data_ == 0) {
		this->data_ = 0;
	}
	else {
		// alloc memory and copy data
		this->data_ = new value_type[in.size_];
		memcpy(this->data_, in.data_, (this->size_)*sizeof(T));
	}

	testInvariant();
}

/// Copy from a View.
///
/// \param in View (source).
///
template<class T> 
template<bool isConst>
Marray<T>::Marray
(
	const View<T, isConst>& in
)
{
	if(!NO_ARG_TEST) {
		in.testInvariant();
	}

	this->dimension_ = in.dimension_;
	this->shape_ = in.shape_;
	this->shapeStrides_ = in.shapeStrides_;
	this->size_ = in.size_;
	this->offset_ = 0;
	this->strides_ = in.shapeStrides_; // !
	this->isSimple_ = true;
	this->coordinateOrder_ = in.coordinateOrder_;

	// alloc memory and copy data
	this->data_ = new value_type[in.size_];
	if(in.isSimple_) {
		memcpy(this->data_, in.data_, (in.size_)*sizeof(T));
	}
	else {
		for(size_t j=0; j<this->size_; ++j)	{
			this->data_[j] = in(j);
		}
	}

	testInvariant();
}

#ifdef HAVE_CPP0X_INITIALIZER_LISTS
/// Constructor.
///
/// \param begin Shape given as initializer list.
/// \param value Value with which all entries are initialized.
/// \param coordinateOrder Flag specifying whether FirstMajorOrder or
/// LastMajorOrder is to be used.
///
template<class T>
Marray<T>::Marray
(
	std::initializer_list<size_t> shape,
	const T& value,
	const CoordinateOrder& coordinateOrder
)
{
	size_t size = 1;
	for(const size_t* it = shape.begin(); it != shape.end(); ++it) {
		Assert(NO_ARG_TEST || static_cast<size_t>(*it) > 0);
		size *= static_cast<size_t>(*it);
	}
	Assert(NO_ARG_TEST || size != 0);
	this->data_ = new value_type[size];
	base::assign(shape.begin(), shape.end(), this->data_, coordinateOrder,
		coordinateOrder); 
	for(size_t j=0; j<size; ++j) {
		this->data_[j] = value;
	}
	testInvariant();
}
#endif

/// Construct Marray with initialization.
///
/// \param begin Iterator to the beginning of a sequence that determines
/// the shape.
/// \param end Iterator to the end of that sequence.
/// \param value Value with which all entries are initialized.
/// \param coordinateOrder Flag specifying whether FirstMajorOrder or
/// LastMajorOrder is to be used.
///
template<class T> 
template<class ShapeIterator>
Marray<T>::Marray
(
	ShapeIterator begin,
	ShapeIterator end,
	const T& value,
	const CoordinateOrder& coordinateOrder
)
{
	size_t size = 1;
	ShapeIterator it = begin;
	while(it != end) {
		Assert(NO_ARG_TEST || static_cast<size_t>(*it) > 0);
		size *= static_cast<size_t>(*it);
		++it;
	}
	Assert(NO_ARG_TEST || size != 0);
	this->data_ = new value_type[size];
	base::assign(begin, end, this->data_, coordinateOrder, coordinateOrder); 
	for(size_t j=0; j<size; ++j) {
		this->data_[j] = value;
	}
	testInvariant();
}

/// Construct Marray without initialization.
///
/// \param is Flag to be set to SkipInitialization.
/// \param begin Iterator to the beginning of a sequence that determines
/// the shape.
/// \param end Iterator to the end of that sequence.
/// \param coordinateOrder Flag specifying whether FirstMajorOrder or
/// LastMajorOrder is to be used.
///
template<class T> 
template<class ShapeIterator>
Marray<T>::Marray
(
	const InitializationSkipping& is,
	ShapeIterator begin,
	ShapeIterator end,
	const CoordinateOrder& coordinateOrder
)
{
	size_t size = 1;
	ShapeIterator it = begin;
	while(it != end) {
		Assert(NO_ARG_TEST || static_cast<size_t>(*it) > 0);
		size *= static_cast<size_t>(*it);
		++it;
	}
	Assert(NO_ARG_TEST || size != 0);
	this->data_ = new value_type[size];
	base::assign(begin, end, this->data_, coordinateOrder, coordinateOrder); 
	testInvariant();
}

/// Destructor.
///
template<class T> 
inline
Marray<T>::~Marray()
{
	delete[] this->data_;
}

/// Assignment.
/// 
/// This operator works as follows:
/// - It always attempts to copy the data from 'in'.
/// - If 'in' and *this have the same size, already allocated memory 
///   is re-used. Otherwise, the memory allocated for *this is freed, 
///   and new memory is allocated to take the copy of 'in'.
/// - If 'in' is un-initialized, memory allocated for *this is freed.
/// - In case of self-assignment, no operation is performed.
/// .
///
/// \param in Marray (source).
/// 
template<class T> 
Marray<T>&
Marray<T>::operator=
(
	const Marray<T>& in
)
{
	testInvariant();
	if(!NO_ARG_TEST) {
		in.testInvariant();
	}
	if(this != &in) { // no self-assignment
		if(in.data_ == 0) {
			delete[] this->data_;
			this->data_ = 0;
			this->dimension_ = 0;
			this->shape_ = coordinate_tuple();
			this->shapeStrides_ = offset_tuple();
			this->size_ = 0;
			this->offset_ = 0;
			this->strides_ = offset_tuple();
			this->isSimple_ = true;
			this->coordinateOrder_ = in.coordinateOrder_;
		}
		else {
			// re-alloc memory if necessary
			if(this->size_ != in.size_) {
				delete[] this->data_;
				this->data_ = new value_type[in.size_];
			}
			this->dimension_ = in.dimension_;
			this->shape_ = in.shape_;
			this->shapeStrides_ = in.shapeStrides_;
			this->size_ = in.size_;
			this->offset_ = 0;
			this->strides_ = in.strides_;
			this->isSimple_ = true;
			this->coordinateOrder_ = in.coordinateOrder_;
			// copy data
			memcpy(this->data_, in.data_, (this->size_)*sizeof(T));
		}
	}
	testInvariant();
	return *this;
}

/// Assignment from View.
///
/// This operator works as follows:
/// - It always attempts to copy the data from 'in'.
/// - If 'in' and *this have overlap, a copy of 'in' is made and 
///   assigned to *this.
/// - If 'in' and *this have the same size, already allocated memory 
///   is re-used. Otherwise, the memory allocated for *this is freed, 
///   and new memory is allocated to take the copy of 'in'.
/// - If 'in' is un-initialized, memory allocated for *this is freed.
/// - In case of self-assignment, no operation is performed.
/// .
/// 
/// \param in View (source).
/// 
template<class T> 
template<bool isConst>
Marray<T>&
Marray<T>::operator=
(
	const View<T, isConst>& in
)
{
	testInvariant();
	if(!NO_ARG_TEST) {
		in.testInvariant();
	}
	if( (void*)(this) != (void*)(&in) ) { // no self-assignment
		if(in.data_ == 0) {
			delete[] this->data_;
			this->data_ = 0;
			this->dimension_ = 0;
			this->shape_ = coordinate_tuple();
			this->shapeStrides_ = offset_tuple();
			this->size_ = 0;
			this->offset_ = 0;
			this->strides_ = offset_tuple();
			this->isSimple_ = true;
			this->coordinateOrder_ = defaultOrder;
		}
		else {
			if(this->overlaps(in)) {
				Marray<T> m = in; // temporary copy
				(*this) = m;
			}
			else {
				// re-alloc memory if necessary
				if(this->size_ != in.size_) {
					delete[] this->data_;
					this->data_ = new value_type[in.size_];
				}
				this->dimension_ = in.dimension_;
				this->shape_ = in.shape_;
				this->shapeStrides_ = in.shapeStrides_;
				this->size_ = in.size_;
				this->offset_ = 0;
				this->strides_ = in.shapeStrides_; // !
				this->isSimple_ = true;
				this->coordinateOrder_ = in.coordinateOrder_;
				if(in.isSimple()) {
					memcpy(this->data_, in.data_, (in.size_)*sizeof(T));
				}
				else {
					for(size_t j=0; j<this->size_; ++j) {
						this->data_[j] = in(j);
					}
				}
			}
		}
	}
	testInvariant();
	return *this;
}

#ifdef HAVE_CPP0X_INITIALIZER_LISTS
/// Resize (existing entries are preserved).
///
/// \param shape Shape given as initializer list.
/// \param value Initial value to be assigned to newly allocated entries.
///
template<class T> 
inline void
Marray<T>::resize
(
	std::initializer_list<size_t> shape,
	const T& value
)
{
	resize(shape.begin(), shape.end(), value);
}
#endif

/// Resize (existing entries are preserved).
///
/// \param begin Iterator to the beginning of a sequence that determines
/// the new shape.
/// \param end Iterator to the end of that sequence.
/// \param value Initial value to be assigned to newly allocated entries.
///
template<class T> 
template<class ShapeIterator>
void
Marray<T>::resize
(
	ShapeIterator begin,
	ShapeIterator end,
	const T& value
)
{	
	testInvariant();

	// compute size
	ShapeIterator it = begin;
	std::vector<size_t> newShape;
	size_t newSize = 1;
	while(it != end) {
		size_t x = static_cast<size_t>(*it);
		Assert(NO_ARG_TEST || x > 0);
		newShape.push_back(x);
		newSize *= static_cast<size_t>(x);
		++it;
	}
	
	// allocate new
	value_type* newData = new value_type[newSize]; 

	// fill all entries with new value
	for(size_t j=0; j<newSize; ++j) {
		newData[j] = value;
	}

	// copy old data in region of overlap
	if(this->data_ != 0) {
		if(newSize == 1) {
			newData[0] = this->data_[0];
		}
		else {
			if(this->dimension_ == 0) {
				newData[0] = this->data_[0];
			}
			else { // neither the new, nor the old marray are scalars
				std::vector<size_t> base1(this->dimension_);
				std::vector<size_t> base2(newShape.size());
				std::vector<size_t> shape1(this->dimension_, static_cast<size_t>(1));
				std::vector<size_t> shape2(newShape.size(), static_cast<size_t>(1));
				for(size_t j=0; j<std::min(static_cast<size_t>(this->dimension_), newShape.size()); ++j) {
					shape1[j] = std::min(this->shape_[j], newShape[j]);
					shape2[j] = shape1[j];
				}
				View<T, true> view1;
				this->constView(base1.begin(), shape1.begin(), view1);
				View<T, false> viewT(newShape.begin(), newShape.end(),
					newData, this->coordinateOrder_,
					this->coordinateOrder_);
				View<T, false> view2;
				viewT.view(base2.begin(), shape2.begin(), view2);
				view1.squeeze();
				view2.squeeze();
				view2 = view1; // copy
			}
		}
		delete[] this->data_; 
		this->data_ = 0;
	}

	base::assign(begin, end, newData, this->coordinateOrder_,
		this->coordinateOrder_);

	testInvariant();
}

/// Invariant test.
///
template<class T> 
inline void
Marray<T>::testInvariant() const
{
	View<T, false>::testInvariant();
	Assert(NO_DEBUG || this->isSimple_);
}

// iterator implementation

/// Empty constructor.
///
/// Sets the view pointer to 0.
///
template<class T, bool isConst>
inline Iterator<T, isConst>::Iterator()
: view_(0), index_(0)
{
}

/// Construct from View on constant data.
///
/// \param view View
/// \param index Index into the View.
///
/// Note for developers: If isConst==false, the construction view_(&view)
/// fails due to incompatible types. This is intended because it should 
/// not be possible to construct a mutable iterator on constant data.
///
template<class T, bool isConst>
inline Iterator<T, isConst>::Iterator
(
	const View<T, true>& view,
	const size_t& index
)
: view_(&view), index_(index)
{
}

/// Construct from View on mutable data.
///
/// \param view View
/// \param index Index into the View.
///
/// Note for developers: If isConst==true, the construction
/// view_(reinterpret_cast<view_pointer>(&view)) works as well.
/// This is intended because it should be possible to construct 
/// a constant iterator on mutable data.
///
template<class T, bool isConst>
inline Iterator<T, isConst>::Iterator
(
	const View<T, false>& view,
	const size_t& index
)
: view_(reinterpret_cast<view_pointer>(&view)), index_(index)
{
}

/// Construct from View on mutable data.
///
/// \param view View
/// \param index Index into the View.
///
/// Note for developers: If isConst==true, the construction
/// view_(reinterpret_cast<view_pointer>(&view)) works as well.
/// This is intended because it should be possible to construct 
/// a constant iterator on mutable data.
///
template<class T, bool isConst>
inline Iterator<T, isConst>::Iterator
(
	View<T, false>& view,
	const size_t& index
)
: view_(reinterpret_cast<view_pointer>(&view)), index_(index)
{
}

/// Copy constructor or conversion from an Iterator on mutable data.
/// 
template<class T, bool isConst>
inline Iterator<T, isConst>::Iterator
(
	const Iterator<T, false>& in
)
: view_(view_pointer(in.view_)), index_(in.index_)
{
}

/// De-reference.
///
template<class T, bool isConst>
inline typename Iterator<T, isConst>::reference
Iterator<T, isConst>::operator*() const
{
	Assert(NO_DEBUG || (view_ != 0 && index_ < view_->size()));
	return (*view_)(index_);
}

/// Pointer.
///
template<class T, bool isConst>
inline typename Iterator<T, isConst>::pointer
Iterator<T, isConst>::operator->() const
{
	Assert(NO_DEBUG || (view_ != 0 && index_ < view_->size()));
	return &((*view_)(index_));
}

/// Element access.
///
template<class T, bool isConst>
inline typename Iterator<T, isConst>::reference
Iterator<T, isConst>::operator[]
(
	const size_t& x
) const
{
	Assert(NO_DEBUG || (view_ != 0 && x+index_ < view_->size()));
	return (*view_)(x+index_);
}

template<class T, bool isConst>
inline Iterator<T, isConst>&
Iterator<T, isConst>::operator+=
(
	const difference_type& x
)
{
	Assert(NO_DEBUG || view_ != 0);
	index_ += x;
	return *this;
}

template<class T, bool isConst>
inline Iterator<T, isConst>&
Iterator<T, isConst>::operator-=
(
	const difference_type& x
)
{
	Assert(NO_DEBUG || view_ != 0);
	index_ -= x;
	return *this;
}

/// Prefix increment.
///
template<class T, bool isConst>
inline Iterator<T, isConst>&
Iterator<T, isConst>::operator++()
{
	Assert(NO_DEBUG || view_ != 0);
	++index_;
	return *this;
}

/// Prefix decrement.
///
template<class T, bool isConst>
inline Iterator<T, isConst>&
Iterator<T, isConst>::operator--()
{
	Assert(NO_DEBUG || view_ != 0);
	Assert(NO_ARG_TEST || index_ > 0);
	--index_;
	return *this;
}

/// Postfix increment.
///
template<class T, bool isConst>
inline Iterator<T, isConst> 
Iterator<T, isConst>::operator++(int)
{
	Assert(NO_DEBUG || view_ != 0);
	Iterator<T, isConst> copy = *this;
	++index_;
	return copy;
}

/// Postfix decrement.
///
template<class T, bool isConst>
inline Iterator<T, isConst> 
Iterator<T, isConst>::operator--(int)
{
	Assert(NO_DEBUG || (view_ != 0 && index_ > 0));
	Iterator<T, isConst> copy = *this;
	--index_;
	return copy;
}

template<class T, bool isConst>
inline Iterator<T, isConst>
Iterator<T, isConst>::operator+
(
	const difference_type& x
) const
{
	Assert(NO_DEBUG || view_ != 0);
	Iterator<T, isConst> tmp = *this;
	tmp += x;
	return tmp;
}

template<class T, bool isConst>
inline Iterator<T, isConst>
Iterator<T, isConst>::operator-
(
	const difference_type& x
) const
{
	Assert(NO_DEBUG || view_ != 0);
	Iterator<T, isConst> tmp = *this;
	tmp -= x;
	return tmp;
}

template<class T, bool isConst>
template<bool isConstLocal>
inline typename Iterator<T, isConst>::difference_type
Iterator<T, isConst>::operator-
(
	const Iterator<T, isConstLocal>& it
) const
{
	Assert(NO_DEBUG || view_ != 0);
	Assert(NO_ARG_TEST || it.view_ != 0);
	return difference_type(index_)-difference_type(it.index_);
}

template<class T, bool isConst>
template<bool isConstLocal> 
inline bool
Iterator<T, isConst>::operator==
(
	const Iterator<T, isConstLocal>& it
) const
{
	Assert(NO_DEBUG || view_ != 0);
	Assert(NO_ARG_TEST || (it.view_ != 0 && (void*)it.view_ == (void*)view_));
	return index_ == it.index_;
}

template<class T, bool isConst>
template<bool isConstLocal>
inline bool
Iterator<T, isConst>::operator!=
(
	const Iterator<T, isConstLocal>& it
) const
{
	Assert(NO_DEBUG || view_ != 0);
	Assert(NO_ARG_TEST || (it.view_ != 0 && it.view_ == view_));
	return index_ != it.index_;
}

template<class T, bool isConst>
template<bool isConstLocal>
inline bool
Iterator<T, isConst>::operator<
(
	const Iterator<T, isConstLocal>& it
) const
{
	Assert(NO_DEBUG || view_ != 0);
	Assert(NO_ARG_TEST || (it.view_ != 0 && it.view_ == view_));
	return(index_ < it.index_);
}

template<class T, bool isConst>
template<bool isConstLocal>
inline bool
Iterator<T, isConst>::operator>
(
	const Iterator<T, isConstLocal>& it
) const
{
	Assert(NO_DEBUG || view_ != 0);
	Assert(NO_ARG_TEST || (it.view_ != 0 && it.view_ == view_));
	return(index_ > it.index_);
}

template<class T, bool isConst>
template<bool isConstLocal>
inline bool
Iterator<T, isConst>::operator<=
(
	const Iterator<T, isConstLocal>& it
) const
{
	Assert(NO_DEBUG || view_ != 0);
	Assert(NO_ARG_TEST || (it.view_ != 0 && it.view_ == view_));
	return(index_ <= it.index_);
}

template<class T, bool isConst>
template<bool isConstLocal> 
inline bool
Iterator<T, isConst>::operator>=
(
	const Iterator<T, isConstLocal>& it
) const
{
	Assert(NO_DEBUG || view_ != 0);
	Assert(NO_ARG_TEST || (it.view_ != 0 && it.view_ == view_));
	return(index_ >= it.index_);
}

/// Fast alternative to comparing with the end iterator.
///
/// \return Boolean indicator.
///
template<class T, bool isConst>
inline bool
Iterator<T, isConst>::hasMore() const
{
	Assert(NO_DEBUG || view_ != 0);
	return index_ < view_->size();
}

/// Get the corresponding index in the View.
///
/// \return index Index.
///
template<class T, bool isConst>
inline size_t 
Iterator<T, isConst>::index() const
{
	return index_;
}

/// Get the corresponding coordinate sequence in the View.
///
/// \param it Iterator into a container starting from which the
/// coordinate sequence is to be written (output).
///
template<class T, bool isConst>
template<class CoordinateIterator>
inline void
Iterator<T, isConst>::coordinate
(
	CoordinateIterator it
) const
{
	Assert(NO_DEBUG || view_ != 0);
	Assert(NO_ARG_TEST || index_ < view_->size());
	view_->indexToCoordinates(index_, it);
}

// Vector implementation

/// Empty constructor.
///
template<class T>
inline
Vector<T>::Vector()
{
	testInvariant();
}

/// Copy constructor.
///
/// \param in Vector (source).
///
template<class T>
template<bool isConst>
Vector<T>::Vector
(
	const View<T, isConst>& in
)
{
	in.testInvariant();
  Assert(NO_ARG_TEST ||
    in.data_ == 0 // un-initialized
    || (in.dimension_ == 0 && in.size_ == 1) // scalar
    || in.dimension_ == 1 // vector
  );
	this->offset_ = 0;
	this->isSimple_ = true;
	this->size_ = in.size();
	this->coordinateOrder_ = in.coordinateOrder();
	if(in.data_ == 0) { // in is uninitialized
		this->data_ = 0;
		this->dimension_ = 0;
		this->shape_ = coordinate_tuple();
		this->shapeStrides_ = offset_tuple();
		this->strides_ = offset_tuple();
	}
	else {
		this->dimension_ = 1;
		this->shape_ = coordinate_tuple(1);
		this->shape_[0] = in.size();
		this->shapeStrides_ = offset_tuple(1);
		this->shapeStrides_[0] = 1;
		this->strides_ = offset_tuple(1);
		this->strides_[0] = 1;
		this->data_ = new T[this->size()];
		if(in.dimension_ == 0) { // in is a scalar
			this->data_[0] = in(0);
		}
		else {
			if(in.isSimple_) {
				memcpy(this->data_, in.data_, (this->size_)*sizeof(T));
			}
			else {
				for(size_t j=0; j<in.size(); ++j) {
					this->data_[j] = in(j);
				}
			}
		}
	}
	testInvariant();
}

/// Construct Vector with initialization.
///
/// \param size Size.
/// \param value Initial value of entries.
/// 
template<class T>
Vector<T>::Vector
(
	const size_t& size,
	const T& value
)
{
  if(size != 0) {
	  size_t shape[1] = {size};
	  this->data_ = new value_type[size];
    base::base::assign(&shape[0], &shape[1], this->data_); 
	  for(size_t j=0; j<size; ++j) {
		  this->data_[j] = value;
	  }
  }
	testInvariant();
}

/// Construct Vector without initialization.
///
/// \param is Flag to be set to SkipInitialization.
/// \param size Size.
/// 
template<class T>
Vector<T>::Vector
(
	const InitializationSkipping& is,
	const size_t& size
) 
{
  if(size != 0) {
	  size_t shape[1] = {size};
	  this->data_ = new value_type[size];
    base::base::assign(&shape[0], &shape[1], this->data_); 
  }
	testInvariant();
}

#ifdef HAVE_CPP0X_INITIALIZER_LISTS
/// Construct with initializer list (C++0x)
///
/// \param list Initializer list.
///
template<class T>
Vector<T>::Vector
(
	std::initializer_list<T> list
)
{
  if(list.size() != 0) {
	  size_t shape[1] = {list.size()};
	  this->data_ = new value_type[list.size()];
    base::base::assign(&shape[0], &shape[1], this->data_);
	  int j=0;
	  for(const T *p = list.begin(); p != list.end(); ++p, ++j) {
		  this->data_[j] = *p;
	  }
  }
	testInvariant();
}
#endif

/// Assignment.
///
/// \param in Vector
///
/// The entries of 'in' are copied.
///
template<class T>
inline Vector<T>&
Vector<T>::operator=
(
	const Vector<T>& in
)
{
	in.testInvariant();
  base::operator=(in);
	testInvariant();
	return *this;
}

/// Assignment from View.
///
/// \param in View.
///
/// The entries of 'in' are copied in the coordinate order of 'in'.
///
template<class T>
template<bool isConst>
inline Vector<T>&
Vector<T>::operator=
(
	const View<T, isConst>& in
)
{
  in.testInvariant();
	Assert(NO_ARG_TEST || 
    in.data_ == 0 || // empty
    (in.dimension_ == 0 && in.size_ == 1) || // scalar
    in.dimension_ == 1); // vector

  if(in.dimension_ == 0 && in.size_ == 1) {
    // in is a View to a scalar
		this->dimension_ = 1;
		this->shape_ = coordinate_tuple(1);
    this->shape_[0] = 1;
		this->shapeStrides_ = offset_tuple(1);
    this->shapeStrides_[0] = 1;
		this->size_ = 1;
		this->offset_ = 0;
		this->strides_ = offset_tuple(1);
    this->strides_[0] = 1;
		this->isSimple_ = true;
		this->coordinateOrder_ = in.coordinateOrder_;
		if(this->size_ != 1) {
			delete[] this->data_;
			this->data_ = new value_type[1];
		}
    this->data_[0] = in(0);
  }
  else {
	  base::operator=(in);
  }
  testInvariant();
  return *this;
}

/// Reshape.
///
/// A Vector can only be reshaped into a Vector which has the same
/// size and thus the same shape. The main pupose of this function
/// is to hide the function reshape inherited from View in a
/// way that is C++ standard complient.
///
/// \param size New size (which has to equal the current size).
///
template<class T>
inline void
Vector<T>::reshape
(
	const size_t& size
)
{
	Assert(size == this->size_);
}

/// Resize.
///
/// \param size New size.
/// \param value Value to be assigned to newly allocated entries.
///
template<class T>
inline void
Vector<T>::resize
(
	const size_t& size,
	const T& value
)
{
  if(size == 0) {
    base::assign();
  }
  else {
	  size_t shape[1] = {size};
	  base::resize(&shape[0], &shape[1], value); 
  }
  testInvariant();
}

/// Element access.
///
template<class T>
inline T&
Vector<T>::operator[]
(
	const size_t& index
)
{
	testInvariant();
	return this->operator()(index);
}

/// Element access.
///
template<class T>
inline const T&
Vector<T>::operator[]
(
	const size_t& index
) const
{
	testInvariant();
	return this->operator()(index);
}

/// Invariant test.
///
template<class T> 
inline void
Vector<T>::testInvariant() const
{
	View<T, false>::testInvariant();
	Assert(NO_DEBUG || 
    this->data_ == 0 ||
		(this->offset_ == 0 && this->isSimple_ && this->dimension_ == 1) 
  );
}

// Matrix implementation

/// Empty constructor.
///
template<class T>
inline
Matrix<T>::Matrix()
{
	testInvariant();
}

/// Copy from a View.
///
/// \param in View (source).
///
template<class T>
template<bool isConst>
Matrix<T>::Matrix
(
	const View<T, isConst>& in
)
{
	in.testInvariant();
	Assert(NO_ARG_TEST || 
    in.data_ == 0 || // not initialized
    (in.dimension_ == 0 && in.size_ == 1) || // scalar
    in.dimension_ == 2); // matrix

	this->offset_ = 0;
	this->isSimple_ = true;
	this->size_ = in.size();
	this->coordinateOrder_ = in.coordinateOrder_;
	if(in.data_ == 0) { // 'in' is uninitialized
		this->data_ = 0;
		this->dimension_ = 0;
		this->shape_ = coordinate_tuple();
		this->shapeStrides_ = offset_tuple();
		this->strides_ = offset_tuple();
	}
	else {
		this->dimension_ = 2;
		this->shape_ = coordinate_tuple(2);
		this->shapeStrides_ = offset_tuple(2);
		this->strides_ = offset_tuple(2);
		if(in.dimension() == 0) { // in is a scalar
			this->shape_[0] = 1;
			this->shape_[1] = 1;
			this->shapeStrides_[0] = 1;
			this->shapeStrides_[1] = 1;
			this->strides_[0] = 1;
			this->strides_[1] = 1;
			this->data_ = new T[1];
			this->data_[0] = in(0);
		}
		else {
			this->shape_[0] = in.shape_[0];
			this->shape_[1] = in.shape_[1];
			this->shapeStrides_[0] = in.shapeStrides_[0];
			this->shapeStrides_[1] = in.shapeStrides_[1];
			this->strides_[0] = in.shapeStrides_[0]; // !
			this->strides_[1] = in.shapeStrides_[1]; // !
			this->data_ = new T[this->size_];
			if(in.isSimple_) {
				memcpy(this->data_, in.data_, (this->size_)*sizeof(T));
			}
			else {
				for(size_t j=0; j<in.size(); ++j) {
					this->data_[j] = in(j);
				}
			}
		}
	}
	testInvariant();
}

/// Construct Matrix with initialization.
///
/// \param n1 size in 1st dimension.
/// \param n2 size in 2nd dimension.
/// \param value Initial value of entries.
/// \param coordinateOrder Flag specifying whether FirstMajorOrder or
/// LastMajorOrder is to be used.
///
template<class T>
inline
Matrix<T>::Matrix
(
	const size_t& n1,
	const size_t& n2,
	const T& value,
	const CoordinateOrder& coordinateOrder
)
{
  if(n1 > 0 && n2 > 0) {
	  size_t shape[2] = {n1, n2};
	  this->data_ = new value_type[n1*n2];
    base::base::assign(&shape[0], &shape[2], this->data_,
      coordinateOrder, coordinateOrder);
	  for(size_t j=0; j<this->size(); ++j) {
		  this->data_[j] = value;
	  }
  }
	testInvariant();
}

/// Construct Matrix without initialization.
///
/// \param is Flag to be set to SkipInitialization.
/// \param n1 Size in 1st dimension.
/// \param n2 Size in 2nd dimension.
/// \param coordinateOrder Flag specifying whether FirstMajorOrder or
/// LastMajorOrder is to be used.
///
template<class T>
inline
Matrix<T>::Matrix
(
	const InitializationSkipping& is, 
	const size_t& n1,
	const size_t& n2,
	const CoordinateOrder& coordinateOrder
) 
{
  if(n1 > 0 && n2 > 0) {
	  size_t shape[2] = {n1, n2};
	  this->data_ = new value_type[n1*n2];
    base::base::assign(&shape[0], &shape[2], this->data_, 
      coordinateOrder, coordinateOrder);
  }
	testInvariant();
}

/// Assignment.
///
/// \param in Matrix (source).
///
template<class T>
inline Matrix<T>&
Matrix<T>::operator=
(
	const Matrix<T>& in
)
{
  in.testInvariant();
	base::operator=(in);
	testInvariant();
	return *this;
}

/// Assignment from View
///
/// \param in View (source).
///
template<class T>
template<bool isConst>
inline Matrix<T>&
Matrix<T>::operator=
(
	const View<T, isConst>& in
)
{
	Assert(NO_ARG_TEST || 
    in.data_ != 0 || // empty
    (in.dimension() == 0 && in.size_ == 1) || // scalar
    in.dimension() == 2);

  if(in.dimension_ == 0 && in.size_ == 1) {
    // in is a View to a scalar
		this->dimension_ = 2;
		this->shape_ = coordinate_tuple(2);
    this->shape_[0] = 1;
    this->shape_[1] = 1;
		this->shapeStrides_ = offset_tuple(2);
    this->shapeStrides_[0] = 1;
    this->shapeStrides_[1] = 1;
		this->size_ = 1;
		this->offset_ = 0;
		this->strides_ = offset_tuple(2);
    this->strides_[0] = 1;
    this->strides_[1] = 1;
		this->isSimple_ = true;
		this->coordinateOrder_ = in.coordinateOrder_;
		if(this->size_ != 1) {
			delete[] this->data_;
			this->data_ = new value_type[1];
		}
    this->data_[0] = in(0);
  }
  else {
	  base::operator=(in);
  }
  testInvariant();
  return *this;
}

/// Resize.
///
/// \param n1 New extension in 1st dimension.
/// \param n2 New extension in 1st dimension.
/// \param value Initial value for newly allocated entries.
///
template<class T>
inline void
Matrix<T>::resize
(
	const size_t& n1,
	const size_t& n2,
	const T& value
)
{
  if(n1 == 0 || n2 == 0) {
    base::assign();
  }
  else {
	  size_t shape[2] = {n1, n2};
	  base::resize(&shape[0], &shape[2], value);
  }
  testInvariant();
}

/// Reshape.
///
/// \param n1 New extension in 1st dimension.
/// \param n2 New extension in 1st dimension.
///
/// Is is necessary that n1*n2 == size().
///
template<class T>
inline void
Matrix<T>::reshape
(
	const size_t& n1,
	const size_t& n2
)
{
	Assert(NO_ARG_TEST || (n2 > 0 && n1 > 0));
	size_t shape[2] = {n1, n2};
	base::reshape(&shape[0], &shape[2]); 
	testInvariant();
}

/// Invariant test.
///
template<class T> 
inline void
Matrix<T>::testInvariant() const
{
	View<T, false>::testInvariant();
	Assert(NO_DEBUG || this->data_ == 0 ||
		(this->offset_ == 0 && this->isSimple_ && this->dimension_ == 2));
}

} // namespace marray 

#endif // #ifndef MARRAY_HXX
