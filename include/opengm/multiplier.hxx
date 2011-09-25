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
#ifndef OPENGM_MULTIPLIER_HXX
#define OPENGM_MULTIPLIER_HXX

namespace opengm {

struct Multiplier {
    // operation (in-place)
    template<class T1, class T2>
    static void op(const T1& in, T2& out)  
    { out *= in; }

    // operation (not in-place)
    template<class T1,class T2,class T3>
    static void op(const T1& in1, const T2& in2, T3& out)  
    { out = in1 * in2; }

    // inverse operation (in-place)
    template<class T1, class T2>
    static void iop(const T1& in, T2& out)  
    { out /= in; }

    // inverse operation (not in-place)
    template<class T1, class T2, class T3>
    static void iop(const T1& in1, const T2& in2, T3& out)  
    { out = in1 / in2; }

    // hyper-operation (in-place)
    template<class T1, class T2>
    static void hop(const T1& in, T2& out) 
    { out.pow(in);  }

    // hyper-operation (not in-place)
    template<class T1, class T2, class T3>
    static void hop(const T1& in1, const T2& in2, T3& out)  
    { out = pow(in1,in2); }

    // inverse hyper-operation (in-place)
    template<class T1,class T2>
    static void ihop(const T1& in, T2& out)  
    { out.pow(1/in); }

    // inverse hyper-operation (not in-place)
    template<class T1, class T2, class T3>
    static void ihop(const T1& in1, const T2& in2, T3& out)  
    { out = pow(in1,1/in2); }

    // neutral element (with return)
    template<class T>
    static T neutral()
    { return static_cast<T>(1); }

    // neutral element (call by reference)
    template<class T>
    static void neutral(T& out)
    { out = static_cast<T>(1); }

    // normalize factor
    template<class T>
    static void normalize(T& out)
    { out.normalize(); }

    // weighted mean
    template<class T1, class T2>
    static void weightedMean(const T1& in1, const T1& in2, const T2& w, T1& out)
    {
        assert(&out != &in2);
        out = in1;
        out /= in2;
        out.pow(w);
        out *= in2;
        // in1^w * in2^(1-w)= (in1/in2)^w*in2 
    }

    // weighted operation
    template<class T1, class T2>
    static void wop(const T2& in, const T1& w, T2& out)
    {
        out *= pow(in,w);
    }

    // inverse weighted operation
    template<class T1, class T2>
    static void iwop(const T2& in, const T1& w, T2& out)
    {
        out *= pow(in,1/w);
    }

    // weighted inverse operation
    template<class T1, class T2>
    static void wiop(const T2& in, const T1& w, T2& out)
    {
        out /= pow(in,w);
    }

    // inverse weighted inverse operation
    template<class T1, class T2>
    static void iwiop(const T2& in, const T1& w, T2& out)
    {
        out /= pow(in,1/w);
    }
};

} // namespace opengm

#endif // #ifndef OPENGM_MULTIPLIER_HXX
