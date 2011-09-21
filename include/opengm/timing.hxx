/// OpenGM. Copyright (c) 2010 by Bjoern Andres and Joerg Hendrik Kappes.
///
/// This software was developed by Bjoern Andres and Joerg Hendrik Kappes.
/// Enquiries shall be directed to:
/// bjoern.andres@iwr.uni-heidelberg.de, kappes@math.uni-heidelberg.de
///
/// Author(s) of this file: Ullrich Koethe
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

#ifndef VIGRA_TIMING_HXX
#define VIGRA_TIMING_HXX

#ifndef VIGRA_NO_TIMING

#include <iostream>
#include <sstream>
#ifdef MULTI_TICTOC
    #include <vector>
#endif
// usage:
// void time_it()
// {
//     USETICTOC;
//     TIC;
//      ...
//     std::cerr << TOC << " for time_it\n";
// }

#ifdef WIN32

    #include "windows.h"

    namespace {

    inline double queryTimerUnit()
    {
        LARGE_INTEGER frequency;
        QueryPerformanceFrequency(&frequency);
        return 1000.0 / frequency.QuadPart;
    }

    inline double tic_toc_diff_num(LARGE_INTEGER const & tic)
    {
        LARGE_INTEGER toc;
        QueryPerformanceCounter(&toc);
        static double unit = queryTimerUnit();
        return ((toc.QuadPart - tic.QuadPart) * unit);
    }

    inline std::string tic_toc_diff_string(LARGE_INTEGER const & tic)
    {
        double diff = tic_toc_diff_num(tic); 
        std::stringstream s;
        s << diff << " msec";
        return s.str();
    }

    inline void tic_toc_diff(LARGE_INTEGER const & tic)
    {
        std::cerr << tic_toc_diff_string(tic) <<std::endl;
    }

    } // unnamed namespace
    
#ifndef MULTI_TICTOC
    #define USETICTOC LARGE_INTEGER tic_timer
    #define TIC QueryPerformanceCounter(&tic_timer)
    #define TOC  tic_toc_diff       (tic_timer)
    #define TOCN tic_toc_diff_num   (tic_timer)
    #define TOCS tic_toc_diff_string(tic_timer)
#else
    #define USETICTOC std::vector<LARGE_INTEGER> tic_timer
    #define TIC tic_timer.push_back(LARGE_INTEGER());\
                QueryPerformanceCounter(&(tic_timer.back()));
    #define TOC  tic_toc_diff       (tic_timer.back());\
                 tic_timer.pop_back();
    #define TOCN tic_toc_diff_num   (tic_timer.back());\
                 tic_timer.pop_back();
    #define TOCS tic_toc_diff_string(tic_timer.back());\
                 tic_timer.pop_back();
#endif

#else

    #if defined(VIGRA_HIRES_TIMING) && !defined(__CYGWIN__)
        // requires linking against librt
    
        #include <time.h>

        namespace {

        inline double tic_toc_diff_num(timespec const & tic)
        {
            timespec toc;
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &toc);
            return ((toc.tv_sec*1000.0 + toc.tv_nsec/1000000.0) -
                  (tic.tv_sec*1000.0 + tic.tv_nsec/1000000.0));
        }

        inline std::string tic_toc_diff_string(timeval const & tic)
        {
            double diff = tic_toc_diff_num(tic); 
            std::stringstream s;
            s << diff << " msec";
            return s.str();
        }

        inline void tic_toc_diff(timeval const & tic)
        {
            std::cerr << tic_toc_diff_string(tic) << std::endl;
        }
        } // unnamed namespace

#ifndef MULTI_TICTOC
        #define USETICTOC timespec tic_timer
        #define TIC clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tic_timer)
        #define TOC  tic_toc_diff       (tic_timer)
        #define TOCN tic_toc_diff_num   (tic_timer)
        #define TOCS tic_toc_diff_string(tic_timer)
#else

        #define USETICTOC std::vector<timespec> tic_timer;
        #define TIC tic_timer.push_back(timespec());\
                    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(tic_timer.back()));
        #define TOC  tic_toc_diff       (tic_timer.back());\
                     tic_timer.pop_back();
        #define TOCN tic_toc_diff_num   (tic_timer.back());\
                     tic_timer.pop_back();
        #define TOCS tic_toc_diff_string(tic_timer.back());\
                     tic_timer.pop_back();
#endif
    #else
    
        #include <sys/time.h>

        namespace {

        inline double tic_toc_diff_num(timeval const & tic)
        {
            timeval toc;
            gettimeofday(&toc, NULL);
            return  ((toc.tv_sec*1000.0 + toc.tv_usec/1000.0) -
                        (tic.tv_sec*1000.0 + tic.tv_usec/1000.0));
        }
        inline std::string tic_toc_diff_string(timeval const & tic)
        {
            double diff = tic_toc_diff_num(tic); 
            std::stringstream s;
            s << diff << " msec";
            return s.str();
        }
        inline void tic_toc_diff(timeval const & tic)
        {
            std::cerr << tic_toc_diff_string(tic)<< std::endl;
        }

        } // unnamed namespace

#ifndef MULTI_TICTOC
        #define USETICTOC timeval tic_timer
        #define TIC  gettimeofday       (&tic_timer, NULL)
        #define TOC  tic_toc_diff       (tic_timer)
        #define TOCN tic_toc_diff_num   (tic_timer)
        #define TOCS tic_toc_diff_string(tic_timer)
#else

        #define USETICTOC std::vector<timeval> tic_timer;
        #define TIC tic_timer.push_back(timeval());\
                    gettimeofday(&(tic_timer.back()), NULL);
        #define TOC  tic_toc_diff       (tic_timer.back());\
                     tic_timer.pop_back();
        #define TOCN tic_toc_diff_num   (tic_timer.back());\
                     tic_timer.pop_back();
        #define TOCS tic_toc_diff_string(tic_timer.back());\
                     tic_timer.pop_back();
#endif

    #endif // VIGRA_HIRES_TIMING

#endif // WIN32




#else // NDEBUG

#define USETICTOC 
#define TIC
#define TOC
#define TOCN
#define TICS
#endif // NDEBUG

#endif // VIGRA_TIMING_HXX
