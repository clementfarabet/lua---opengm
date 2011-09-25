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
#ifndef OPENGM_EXITCONDITION_HXX
#define OPENGM_EXITCONDITION_HXX

#include <limits>
#include <time.h>

namespace opengm {

template<class T, class OP, class ACC>
class ExitConditions {
public:
    static const size_t EXITCOND_FLAG_TIME     = 1; // less than t seconds
    static const size_t EXITCOND_FLAG_QUANTITY = 2; // maximal the N best
    static const size_t EXITCOND_FLAG_QUALITY  = 4; // better than alpha
    static const size_t EXITCOND_FLAG_RELBEST  = 8; // not wrose than epsilon to the optimum

    size_t conditionflag;
    float max_timediff;
    size_t max_numberOfSolutions;
    T rel_valueThreshold;
    T abs_valueThreshold;

    ExitConditions();
    void setStartExitCondTime()
        { starttime = clock(); };
    size_t checkExitCond(size_t cond, T value, size_t num_sol);

private:
    clock_t starttime;
    T best_value;
    bool best_value_isSet;
};

template<class T, class OP, class ACC >
ExitConditions<T,OP,ACC>::ExitConditions()
{
    conditionflag = EXITCOND_FLAG_QUANTITY;
    max_timediff = std::numeric_limits<float>::infinity();
    starttime = clock();
    max_numberOfSolutions = 1;
    best_value_isSet = false; 
    ACC::neutral(best_value);
    ACC::neutral(abs_valueThreshold);
    OP::neutral(rel_valueThreshold); 
}

template<class T, class OP, class ACC >
size_t 
ExitConditions<T,OP,ACC>::checkExitCond
(
    size_t cond, 
    T value, 
    size_t num_sol
)
{
    unsigned int ret = 0;
    if(cond & conditionflag & EXITCOND_FLAG_TIME){
        if((float)(clock()-starttime)/CLOCKS_PER_SEC > max_timediff){
            ret = ret |  EXITCOND_FLAG_TIME;
        }
    }
    if(cond & conditionflag & EXITCOND_FLAG_QUANTITY){
        if(num_sol == max_numberOfSolutions){
            ret = ret |  EXITCOND_FLAG_QUANTITY; 
        }
    }
    if(cond & conditionflag & EXITCOND_FLAG_QUALITY){
        if(ACC::ibop(value,abs_valueThreshold)){
            ret = ret | EXITCOND_FLAG_QUALITY;
        } 
    }
    if(cond & conditionflag & EXITCOND_FLAG_RELBEST){
        if(num_sol == 0){
            best_value       = value; 
            best_value_isSet = true;
        }
        else{
            T bound;
            OP::op(best_value,rel_valueThreshold,bound);
            if(ACC::ibop(value,bound)){ //example min-sum : value > best + threshold
                ret = ret | EXITCOND_FLAG_RELBEST;
            }
        }
    }
    return ret;
}

} // namespace opengm

#endif // #ifndef OPENGM_EXITCONDITION_HXX
