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
#ifndef OPENGM_HXX
#define OPENGM_HXX

#include <stddef.h>
#include <stdexcept>
#include <functional>

namespace opengm {

// assertion testing

#ifdef NDEBUG
	const bool NO_DEBUG = true; 
	const bool NO_ARG_TEST = true;
#else
	const bool NO_DEBUG = false;
	const bool NO_ARG_TEST = false;
#endif

template<class A> inline void Assert(A assertion) {
	if(!assertion) throw std::runtime_error("Assertion failed.");
}

// general functors

template<class T>
struct Reciprocal
: public std::unary_function<T, T>
{
	inline T operator()(const T& x)	const
		{ return T(1)/x; }
};

template<class T>
struct AbsoluteDifference
: public std::binary_function<T, T, T>
{
	inline T operator()(const T& x1, const T& x2) const
		{ if(x1 < x2) return x2-x1; else return x1-x2; }
};

} // namespace opengm 

#endif // #ifndef OPENGM_HXX
