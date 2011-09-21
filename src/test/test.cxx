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
#include <vector>
#include <set>
#include <string>
#include <iostream>

#include <opengm/discretespace.hxx>
#include <opengm/explicitfactor.hxx>
#include <opengm/graphicalmodel.hxx>
#include <opengm/adder.hxx>
#include <opengm/accumulation.hxx>
#include <opengm/multiplier.hxx>
#include <opengm/minimizer.hxx>
#include <opengm/integrator.hxx>
#include <opengm/maximizer.hxx>
#include <opengm/maxdistance.hxx>
#include <opengm/serialization.hxx>
#include <opengm/decomposer.hxx>

#include <opengm/inference/beliefpropagation.hxx>
#include <opengm/inference/treereweightedbeliefpropagation.hxx>
#include <opengm/inference/movemaker.hxx>
#include <opengm/inference/astar.hxx>
#include <opengm/inference/icm.hxx>
#include <opengm/inference/lazyflipper.hxx>

#define test(x) \
    if(!(x)) { \
    std::stringstream s; \
    s << #x << " does not hold [line " << __LINE__ << "]"; \
    throw std::logic_error(s.str().c_str()); \
    exit(1); \
}

#define testEqual(x, y) \
    if(!(x == y)) { \
    std::stringstream s; \
    s << x << " != " << y << " [line " << __LINE__ << ": " << #x << " == " << #y << " ]"; \
    throw std::logic_error(s.str().c_str()); \
    exit(1); \
}

#define testEqualTolerance(x, y, epsilon) \
    if( (x<y && y-x > epsilon) || (x>y && x-y > epsilon) ) { \
    std::stringstream s; \
    s << x << " != " << y << " [line " << __LINE__ << ": " << #x << " == " << #y << " ]"; \
    throw std::logic_error(s.str().c_str()); \
    exit(1); \
}

template<class It1, class It2>
void testEqualSequence(It1 begin1, It1 end1, It2 begin2) {
    while(begin1 != end1) {
        test(*begin1 == *begin2);
        ++begin1;
        ++begin2;
    }
}

#define DGM_UNUSED(x) (void)x

struct Evencounter {   
    template<class T>
    static T op(const T& in1, const T& in2) {
        if(static_cast<size_t>(in1) % 2 == 0) 
            return in2 + 1;
        else
            return in2;
    }

    template<class T>
    static bool bop(const T& in1, const T& in2) {
        return false;
    }  

    template<class T>
    static void op(const T& in1, const T& in2, T& out) {
        out = in2; 
        if(static_cast<size_t>(in1) % 2 == 0)
            ++out;
    }

    template<class T>
    static T neutral() {
        return static_cast<T>(0);
    }  

    template<class T>
    static void neutral(T& out) {
        out = static_cast<T>(0);
    }
};

struct ExplicitFactorTest
{
    opengm::DiscreteSpace *space, *otherSpace;
    opengm::ExplicitFactor<float> p0, p1, p2, p3, p4, p5;

    template<class T>
    static void testEqualExplicitFactor
    (
        const opengm::ExplicitFactor<T>& f,
        const opengm::ExplicitFactor<T>& g
    ) 
    {
        test(&f.space() == &g.space());
        test(f.isConstant() == g.isConstant());
        if(!f.isConstant()) {
            test(f.numberOfVariables() == g.numberOfVariables());
            std::vector<size_t> vi(f.numberOfVariables());
            f.variableIndices(vi.begin());
            std::vector<size_t> vi2(g.numberOfVariables());
            g.variableIndices(vi2.begin());
            testEqualSequence(vi.begin(), vi.end(), vi2.begin());
            size_t numberOfEntries = 1;
            for(size_t j=0; j<vi.size(); ++j) {
                numberOfEntries *= f.space().numberOfStates(vi[j]);
            }
            for(size_t j=0; j<numberOfEntries; ++j) {
                test(f(j) == g(j));
            }
        }
        else {
            test(f() == g());
        }
    }

    ExplicitFactorTest() {
        std::set<size_t> variableIndices;
        std::vector<size_t> numbersOfStates;

        // test spaces:
        numbersOfStates.resize(5);
        numbersOfStates[0] = 2;
        numbersOfStates[1] = 3;
        numbersOfStates[2] = 2;
        numbersOfStates[3] = 2;
        numbersOfStates[4] = 1;
        space = new opengm::DiscreteSpace(numbersOfStates.begin(), numbersOfStates.end());
        numbersOfStates.resize(2);
        numbersOfStates[0] = 3;
        numbersOfStates[1] = 2;
        otherSpace = new opengm::DiscreteSpace(numbersOfStates.begin(), numbersOfStates.end());

        // test potentials:
        // p0: 3rd order potential
        variableIndices.insert(0);
        variableIndices.insert(1);
        variableIndices.insert(2);
        p0.assign(*space, variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        p0(0,0,0) = 2;
        p0(1,0,0) = 3;
        p0(0,1,0) = 4;
        p0(1,1,0) = 5;
        p0(0,2,0) = 6;
        p0(1,2,0) = 7;
        p0(0,0,1) = 8;
        p0(1,0,1) = 9;
        p0(0,1,1) = 10;
        p0(1,1,1) = 11;
        p0(0,2,1) = 12;
        p0(1,2,1) = 13;
        // p1: 1st order potential
        variableIndices.insert(1);
        p1.assign(*space, variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        p1(0) = 2;
        p1(1) = 3;
        p1(2) = 4;
        // p2: constant potential
        p2.assign(*space, 2);
        test(p2.isConstant());
        // p3: potential depending on one variable with only one state
        variableIndices.insert(4);
        p3.assign(*space, variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        p3(0) = 2;
        // p4: ternary potential that has a non-trivial variable overlap with p0
        variableIndices.insert(0);
        variableIndices.insert(2);
        variableIndices.insert(3);
        p4.assign(*space, variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        p4(0,0,0) = 2;
        p4(1,0,0) = 3;
        p4(0,1,0) = 4;
        p4(1,1,0) = 5;
        p4(0,0,1) = 6;
        p4(1,0,1) = 7;
        p4(0,1,1) = 8;
        p4(1,1,1) = 9;
        // p5: potential of another discrete space
        variableIndices.insert(0);
        variableIndices.insert(1);
        p5.assign(*otherSpace, variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
    }

    ~ExplicitFactorTest() {
        delete space;
        delete otherSpace;
    }

    void constructionAndQuery(void){
        opengm::ExplicitFactor<float>* p;
        std::set<size_t> variableIndices;
        std::vector<size_t> vi;
        std::vector<size_t> states1(2), states2(3);

        // constant potential constructor:
        p = new opengm::ExplicitFactor<float>(*space, 2);

        test(&(p->space()) == space);
        test(p->isConstant());
        test((*p)() == 2); 
        test(p->numberOfVariables() == 0);
        for(size_t j=0; j<space->dimension()+1; ++j) {
            // note that we test also for invalid variable indices
            std::set<size_t> vi;
            vi.insert(j);
            test(!p->dependsOnVariables(vi.begin(), vi.end()));
            test(!p->dependsOnVariable(j));
        }
        test((*p)(states1.begin()) == (*p)()); // assumed semantics: although the state might be invalid, at() must return the constant!
        delete p;

        // variable indices constructor with non-empty valid set of variable indices:
        variableIndices.insert(1);
        variableIndices.insert(3);
        p = new opengm::ExplicitFactor<float>(*space, variableIndices.begin(), variableIndices.end());

        test(&(p->space()) == space);
        test(!p->isConstant());
        test(p->numberOfVariables() == 2);
        vi.resize(p->numberOfVariables());
        p->variableIndices(vi.begin());
        testEqualSequence(vi.begin(), vi.end(), variableIndices.begin());
        variableIndices.clear();
        variableIndices.insert(1);
        test(p->dependsOnVariables(variableIndices.begin(), variableIndices.end()));
        variableIndices.insert(3);
        test(p->dependsOnVariables(variableIndices.begin(), variableIndices.end()));
        variableIndices.insert(4); // valid
        test(!p->dependsOnVariables(variableIndices.begin(), variableIndices.end()));
        variableIndices.insert(900); // invalid
        test(!p->dependsOnVariables(variableIndices.begin(), variableIndices.end()));
        test(p->dependsOnVariable(1));
        test(p->dependsOnVariable(3));
        test(!p->dependsOnVariable(4)); // valid
        test(!p->dependsOnVariable(900)); // invalid
    }

    void assignment(void)
    {
        opengm::ExplicitFactor<float> p;
        std::set<size_t> variableIndices;
        std::vector<size_t> states(2);
        p = p0; // potential depending on variables

        p.assign(*otherSpace, 2); // make a constant potential

        test(&(p.space()) == otherSpace);
        test(p.isConstant());
        test(p() == 2); 
        test(p.numberOfVariables() == 0);

        p = opengm::ExplicitFactor<float>(*otherSpace, 2); // constant potential
        variableIndices.insert(1);
        variableIndices.insert(3);
        p.assign(*space, variableIndices.begin(), variableIndices.end());

        test(&(p.space()) == space);
        test(!p.isConstant());
        test(p.numberOfVariables() == 2);
        std::vector<size_t> vi(2);
        p.variableIndices(vi.begin());
        test(vi[0] == 1);
        test(vi[1] == 3);
        for(size_t x1 = 0; x1<space->numberOfStates(1); ++x1)
            for(size_t x3 = 0; x3<space->numberOfStates(3); ++x3) {
                states[0] = x1;
                states[1] = x3;
                p(states.begin()) = 2; // test not throw exception
            }
    }

    static float by2plus1(const float& x) {return(2*x+1);}; // a unary operation on floats
    static float multiply(const float& x, const float& y) {return(x*y);} // a binary operation on floats
    void operation(void)
    {
        opengm::ExplicitFactor<float> p;
        std::vector<size_t> expectedVariableIndices(4);
        std::vector<size_t> vi;

        // unary (ExplicitFactor) operation:
        p = p0;
        p.operate(&ExplicitFactorTest::by2plus1);

        test(p(0,0,0) == 1+2*p0(0,0,0));
        test(p(1,0,0) == 1+2*p0(1,0,0));
        test(p(0,1,0) == 1+2*p0(0,1,0));
        test(p(1,1,0) == 1+2*p0(1,1,0));
        test(p(0,2,0) == 1+2*p0(0,2,0));
        test(p(1,2,0) == 1+2*p0(1,2,0));
        test(p(0,0,1) == 1+2*p0(0,0,1));
        test(p(1,0,1) == 1+2*p0(1,0,1));
        test(p(0,1,1) ==1+2*p0(0,1,1));
        test(p(1,1,1) == 1+2*p0(1,1,1));
        test(p(0,2,1) == 1+2*p0(0,2,1));
        test(p(1,2,1) == 1+2*p0(1,2,1));

        p = p1;
        p.operate(&ExplicitFactorTest::by2plus1);

        test(p(0) == 1+2*p1(0));
        test(p(1) == 1+2*p1(1));
        test(p(2) == 1+2*p1(2));

        p = p2;
        p.operate(&ExplicitFactorTest::by2plus1);

        test(p.isConstant());
        test(p() == 1+2*p2());

        p = p3;
        p.operate(&ExplicitFactorTest::by2plus1);

        test(p(0) == 1+2*p3(0));

        // mixed binary (ExplicitFactor, VALUE) operation:

        p = p0;
        p.operate(2, &ExplicitFactorTest::multiply);

        test(p(0,0,0) == 2*p0(0,0,0));
        test(p(1,0,0) == 2*p0(1,0,0));
        test(p(0,1,0) == 2*p0(0,1,0));
        test(p(1,1,0) == 2*p0(1,1,0));
        test(p(0,2,0) == 2*p0(0,2,0));
        test(p(1,2,0) == 2*p0(1,2,0));
        test(p(0,0,1) == 2*p0(0,0,1));
        test(p(1,0,1) == 2*p0(1,0,1));
        test(p(0,1,1) == 2*p0(0,1,1));
        test(p(1,1,1) == 2*p0(1,1,1));
        test(p(0,2,1) == 2*p0(0,2,1));
        test(p(1,2,1) == 2*p0(1,2,1));

        p = p1;
        p.operate(2, &ExplicitFactorTest::multiply);

        test(p(0) == 2*p1(0));
        test(p(1) == 2*p1(1));
        test(p(2) == 2*p1(2));

        p = p2;
        p.operate(2, &ExplicitFactorTest::multiply);

        test(p.isConstant());
        test(p() == 2*p2());

        p = p3;
        p.operate(2, &ExplicitFactorTest::multiply);

        test(p(0) == 2*p3(0));

        // p0 * p4, non-trivial variable overlap
        expectedVariableIndices.resize(4);
        expectedVariableIndices[0] = 0;
        expectedVariableIndices[1] = 1;
        expectedVariableIndices[2] = 2;
        expectedVariableIndices[3] = 3;
        opengm::ExplicitFactor<float>::operate(p0, p4, p, &ExplicitFactorTest::multiply);

        test(!p.isConstant());
        test(p.numberOfVariables() == 4);
        vi.resize(p.numberOfVariables());
        p.variableIndices(vi.begin());
        testEqualSequence(vi.begin(), vi.end(), expectedVariableIndices.begin());
        for(size_t x3=0; x3<2; ++x3)
            for(size_t x2=0; x2<2; ++x2)
                for(size_t x1=0; x1<3; ++x1)
                    for(size_t x0=0; x0<2; ++x0) {
                        test(p(x0,x1,x2,x3) == p0(x0,x1,x2)*p4(x0,x2,x3));
                    }

                    // p0 * p2, non-constant times constant potential:
                    expectedVariableIndices.resize(3);
                    expectedVariableIndices[0] = 0;
                    expectedVariableIndices[1] = 1;
                    expectedVariableIndices[2] = 2;
                    opengm::ExplicitFactor<float>::operate(p0, p2, p, &ExplicitFactorTest::multiply);

                    test(p.numberOfVariables() == 3);
                    vi.resize(p.numberOfVariables());
                    p.variableIndices(vi.begin());
                    testEqualSequence(vi.begin(), vi.end(), expectedVariableIndices.begin());
                    for(size_t x2=0; x2<2; ++x2)
                        for(size_t x1=0; x1<3; ++x1)
                            for(size_t x0=0; x0<2; ++x0)
                            {
                                test(p(x0,x1,x2) == p0(x0,x1,x2)*p2());
                            }

                            // p2 * p0, constant times non-constant potential:
                            expectedVariableIndices.resize(3);
                            expectedVariableIndices[0] = 0;
                            expectedVariableIndices[1] = 1;
                            expectedVariableIndices[2] = 2;
                            opengm::ExplicitFactor<float>::operate(p2, p0, p, &ExplicitFactorTest::multiply);

                            test(p.numberOfVariables() == 3);
                            vi.resize(p.numberOfVariables());
                            p.variableIndices(vi.begin());
                            testEqualSequence(vi.begin(), vi.end(), expectedVariableIndices.begin());
                            for(size_t x2=0; x2<2; ++x2)
                                for(size_t x1=0; x1<3; ++x1)
                                    for(size_t x0=0; x0<2; ++x0) {
                                        test(p(x0,x1,x2) == p0(x0,x1,x2)*p2());
                                    }

                                    // p4 * p1, p1 has 1 variable on which p4 does not depend:
                                    expectedVariableIndices.resize(4);
                                    expectedVariableIndices[0] = 0;
                                    expectedVariableIndices[1] = 1;
                                    expectedVariableIndices[2] = 2;
                                    expectedVariableIndices[3] = 3;
                                    opengm::ExplicitFactor<float>::operate(p4, p1, p, &ExplicitFactorTest::multiply);

                                    test(p.numberOfVariables() == 4);
                                    vi.resize(p.numberOfVariables());
                                    p.variableIndices(vi.begin());
                                    testEqualSequence(vi.begin(), vi.end(), expectedVariableIndices.begin());
                                    for(size_t x3=0; x3<2; ++x3)
                                        for(size_t x2=0; x2<2; ++x2)
                                            for(size_t x1=0; x1<3; ++x1)
                                                for(size_t x0=0; x0<2; ++x0) {
                                                    test(p(x0,x1,x2,x3) == p1(x1)*p4(x0,x2,x3));
                                                }

                                                // p1 * p4, p4 has 3 variables on which p1 does not depend:
                                                expectedVariableIndices.resize(4);
                                                expectedVariableIndices[0] = 0;
                                                expectedVariableIndices[1] = 1;
                                                expectedVariableIndices[2] = 2;
                                                expectedVariableIndices[3] = 3;
                                                opengm::ExplicitFactor<float>::operate(p1, p4, p, &ExplicitFactorTest::multiply);

                                                test(p.numberOfVariables() == 4);
                                                vi.resize(p.numberOfVariables());
                                                p.variableIndices(vi.begin());
                                                testEqualSequence(vi.begin(), vi.end(), expectedVariableIndices.begin());
                                                for(size_t x3=0; x3<2; ++x3)
                                                    for(size_t x2=0; x2<2; ++x2)
                                                        for(size_t x1=0; x1<3; ++x1)
                                                            for(size_t x0=0; x0<2; ++x0) {
                                                                test(p(x0,x1,x2,x3) == p1(x1)*p4(x0,x2,x3));
                                                            }

                                                            // additions. ??? test results!
                                                            // unary +
                                                            p = +p0; 
                                                            // basic types first
                                                            p = 2.0f + p0;
                                                            p = 2.0f - p0;
                                                            p = 2.0f * p0;
                                                            p = 2.0f / p0;
                                                            // exponential operations
                                                            p = p0;
                                                            p.pow(2);

                                                            p = p0;
                                                            p = pow(p, 2);

                                                            p = p0;
                                                            p.pow(2.0f);

                                                            p = p0;
                                                            p = pow(p, 2.0f);

                                                            p = p0;
                                                            p.exp();

                                                            p = p0;
                                                            p = exp(p);

                                                            p = p0;
                                                            p.log();

                                                            p = p0;
                                                            p = log(p);
    }

    void manipulation(void)
    {
        opengm::ExplicitFactor<float> p;
        std::vector<size_t> states(4);
        states[0] = 1;
        states[1] = 2;
        states[2] = 3;
        states[3] = 4;

        std::vector<size_t> assignmentDimensions;
        std::vector<size_t> assignmentStates;

        // assign variables of a constant potential
        p = p2;
        assignmentDimensions.resize(2); // variable and state invalid. assumed semantics: suppress
        assignmentStates.resize(2);
        assignmentDimensions[0] = 1;
        assignmentDimensions[1] = 900;
        assignmentStates[0] = 1;
        assignmentStates[1] = 901;
        p.fixVariables(assignmentDimensions.begin(), assignmentDimensions.end(), assignmentStates.begin());

        test(p.isConstant());
        test(p.numberOfVariables() == 0);
        test(p() == 2);

        // assign one of three variables
        p = p0;
        p.fixVariables(assignmentDimensions.begin(), assignmentDimensions.end(), assignmentStates.begin());

        test(!p.isConstant());
        test(p.numberOfVariables() == 2);
        std::vector<size_t> vi(2);
        p.variableIndices(vi.begin());
        test(vi[0] == 0);
        test(vi[1] == 2);
        test(p(0,0) == 4);
        test(p(1,0) == 5);
        test(p(0,1) == 10);
        test(p(1,1) == 11);

        // assign two of three variables
        p = p0;
        assignmentDimensions.resize(3); 
        assignmentStates.resize(3);
        assignmentDimensions[0] = 0;
        assignmentDimensions[1] = 2;
        assignmentDimensions[2] = 900;
        assignmentStates[0] = 0;
        assignmentStates[1] = 1;
        assignmentStates[2] = 901;
        p.fixVariables(assignmentDimensions.begin(), assignmentDimensions.end(), assignmentStates.begin());

        test(!p.isConstant());
        testEqual(p.numberOfVariables(), 1);
        vi.resize(1);
        p.variableIndices(vi.begin());
        testEqual(vi[0], 1);
        testEqual(p(0),  8);
        testEqual(p(1),  10);
        testEqual(p(2),  12);

        // assign three of three variables
        p = p0;
        assignmentDimensions.resize(4); 
        assignmentStates.resize(4);
        assignmentDimensions[0] = 0;
        assignmentDimensions[1] = 1;
        assignmentDimensions[2] = 2;
        assignmentDimensions[3] = 900;
        assignmentStates[0] = 0;
        assignmentStates[1] = 2;
        assignmentStates[2] = 1;
        assignmentStates[3] = 901;
        p.fixVariables(assignmentDimensions.begin(), assignmentDimensions.end(), assignmentStates.begin());

        test(&(p.space()) == space);
        test(p.isConstant());
        test(p() == 12);
        test(p.numberOfVariables() == 0);
        for(size_t j=0; j<space->dimension()+1; ++j) {
            // note that we test also for invalid variable indices
            std::set<size_t> vi;
            vi.insert(j);
            test(!p.dependsOnVariables(vi.begin(), vi.end()));
            test(!p.dependsOnVariable(j));
        }
        test(p(states.begin()) == p()); // assumed semantics: although the state might be invalid, at() must return the constant!

        // setConstant() test
        p = p0;
        p.assign(p.space(), 2);

        test(&(p.space()) == space);
        test(p.isConstant());
        test(p() == 2); 
        test(p.numberOfVariables() == 0);
        for(size_t j=0; j<space->dimension()+1; ++j) {
            // note that we test also for invalid variable indices
            std::set<size_t> vi;
            vi.insert(j);
            test(!p.dependsOnVariables(vi.begin(), vi.end()));
            test(!p.dependsOnVariable(j));
        }
        test(p(states.begin()) == p()); // assumed semantics: although the state might be invalid, at() must return the constant!
    }

    void accumulation(){      

        //* Testing Full Accumulation
        //****************************
        {
            float v;
            std::vector<size_t> states;

            // potential with three variables
            p0.accumulate<Evencounter>(v,states);
            test(v == 6);
            p0.accumulate<opengm::Maximizer>(v,states);
            test(v == 13);
            p0.accumulate<opengm::Minimizer>(v,states);
            test(v == 2);       

            p1.accumulate<Evencounter>(v,states);
            test(v == 2);

            // constant potential
            p2.accumulate<Evencounter>(v,states);
            test(v == 1);

            // potential with one variable with one state
            p3.accumulate<Evencounter>(v,states);
            test(v == 1);
        }

        //* Testing Partial Accumulation
        //********************************
        {
            std::set<size_t> accVariableIndices; //empty 
            std::set<size_t> variableIndices; //empty

            opengm::ExplicitFactor<float>                values; 
            opengm::ExplicitFactor<std::vector<size_t> > states;

            values = p0;

            // empty set of variable indices,
            p0.accumulate<Evencounter>(accVariableIndices.begin(), 
                accVariableIndices.end(), values);
            testEqual(values(0), 1);
            testEqual(values(1), 0);
            testEqual(values(2), 1);
            testEqual(values(3), 0);
            testEqual(values(4), 1);
            testEqual(values(5), 0);
            testEqual(values(6), 1);
            testEqual(values(7), 0);
            testEqual(values(8), 1);
            testEqual(values(9), 0);
            testEqual(values(10), 1);
            testEqual(values(11), 0);

            // one of three variables
            accVariableIndices.insert(1);
            accVariableIndices.insert(900); // assumed semantics: invalid variable indices are suppressed 
            variableIndices.insert(0);
            variableIndices.insert(2);
            values.assign(*space, variableIndices.begin(), variableIndices.end());

            p0.accumulate<Evencounter>(accVariableIndices.begin(), 
                accVariableIndices.end(), values);

            testEqual(values.numberOfVariables(), 2);
            testEqual(values.variableIndex(0), 0);
            testEqual(values.variableIndex(1), 2);
            testEqual(values(0), 3);
            testEqual(values(1), 0);
            testEqual(values(2), 3);
            testEqual(values(3), 0);

            // two of three variables
            accVariableIndices.clear();   
            accVariableIndices.insert(0);
            accVariableIndices.insert(900); // assumed semantics: invalid variable indices are suppressed
            accVariableIndices.insert(2);
            variableIndices.clear();
            variableIndices.insert(1);    
            values.assign(*space, variableIndices.begin(), variableIndices.end());
            p0.accumulate<Evencounter>(accVariableIndices.begin(), 
                accVariableIndices.end(), values);

            test(values.numberOfVariables() == 1);
            test(values.variableIndex(0) == 1);
            test(values(0) == 2);
            test(values(1) == 2);
            test(values(2) == 2);


            // three of three variables
            accVariableIndices.clear(); 
            accVariableIndices.insert(900); // assumed semantics: invalid variable indices are suppressed
            accVariableIndices.insert(0);
            accVariableIndices.insert(1);
            accVariableIndices.insert(2);
            variableIndices.clear();

            values.assign(*space);
            p0.accumulate<Evencounter>(accVariableIndices.begin(), 
                accVariableIndices.end(), values);

            test(values.numberOfVariables() == 0);
            test(values(0) == 6);


            // constant potential
            accVariableIndices.clear();
            accVariableIndices.insert(900); // assumed semantics: invalid variable indices are suppressed
            accVariableIndices.insert(0);
            accVariableIndices.insert(1); 
            values.assign(*space);
            p2.accumulate<Evencounter>(accVariableIndices.begin(), 
                accVariableIndices.end(), values);

            test(values.numberOfVariables() == 0);
            test(values(0) == 1);
        }
    }

    void maxMinInt() {
        size_t vi[] = {0, 1, 3};
        opengm::ExplicitFactor<float> p(*space, vi, vi+3);

        // maximize, sequence of indices
        {
            size_t vi2[] = {0, 3};
            opengm::ExplicitFactor<float> tmp = p; // copy
            tmp.maximize(vi2, vi2+2);
            // ??? test result
        }
        // maximize, one index
        {
            size_t vi2 = 3;
            opengm::ExplicitFactor<float> tmp = p; // copy
            tmp.maximize(vi2);
            // ??? test result
        }
        // maximize, all variables
        {
            opengm::ExplicitFactor<float> tmp = p; // copy
            tmp.maximize();
            // ??? test result
        }

        // minimize, sequence of indices
        {
            size_t vi2[] = {0, 3};
            opengm::ExplicitFactor<float> tmp = p; // copy
            tmp.minimize(vi2, vi2+2);
            // ??? test result
        }
        // minimize, one index
        {
            size_t vi2 = 3;
            opengm::ExplicitFactor<float> tmp = p; // copy
            tmp.minimize(vi2);
            // ??? test result
        }
        // minimize, all variables
        {
            opengm::ExplicitFactor<float> tmp = p; // copy
            tmp.minimize();
            // ??? test result
        }

        // integrate, sequence of indices
        {
            size_t vi2[] = {0, 3};
            opengm::ExplicitFactor<float> tmp = p; // copy
            tmp.integrate(vi2, vi2+2);
            // ??? test result
        }
        // integrate, one index
        {
            size_t vi2 = 3;
            opengm::ExplicitFactor<float> tmp = p; // copy
            tmp.integrate(vi2);
            // ??? test result
        }
        // integrate, all variables
        {
            opengm::ExplicitFactor<float> tmp = p; // copy
            tmp.integrate();
            // ??? test result
        }
    }

    void normalizeSubtractOffset() {
        {
            size_t vi[] = {0, 1};
            opengm::ExplicitFactor<float> p(*space, vi, vi+2);
            float n = 0;
            for(size_t j=0; j<p.table().size(); ++j) {
                p(j) = static_cast<float>(j);
                n += static_cast<float>(j);
            }
            p.normalize();
            for(size_t j=0; j<p.table().size(); ++j) {
                testEqualTolerance(p(j), static_cast<float>(j)/n, 1e-6);
            }
        }
        {
            size_t vi[] = {0, 1};
            opengm::ExplicitFactor<float> p(*space, vi, vi+2);
            for(size_t j=0; j<p.table().size(); ++j) {
                p(j) = static_cast<float>(j+4);
            }
            p.subtractOffset();
            for(size_t j=0; j<p.table().size(); ++j) {
                testEqual(p(j), static_cast<float>(j));
            }
        }
    }

    void isSubmodularTest() {
        size_t numbersOfStates[] = {2,2};
        opengm::DiscreteSpace space(numbersOfStates, numbersOfStates+2);

        // constant factor
        {
            opengm::ExplicitFactor<float> f(space, 42);
            test(f.isSubmodular());
        }

        // unary factors
        {
            size_t variableIndices[] = {0};
            opengm::ExplicitFactor<float> f(space, variableIndices, variableIndices+1);
            f(0) = 1;
            f(1) = 2;
            test(f.isSubmodular());
            f(0) = 2;
            f(1) = 1;
            test(f.isSubmodular());
            f(0) = 1;
            f(1) = 1;
            test(f.isSubmodular());
        }

        // binary factors
        {
            size_t variableIndices[] = {0,1};
            opengm::ExplicitFactor<float> f(space, variableIndices, variableIndices+2);
            f(0,0) = 1;
            f(0,1) = 1;
            f(1,0) = 1;
            f(1,1) = 1;
            test(f.isSubmodular());
            f(0,0) = 1;
            f(0,1) = 1;
            f(1,0) = 2;
            f(1,1) = 1;
            test(f.isSubmodular());
            f(0,0) = 1;
            f(0,1) = 2;
            f(1,0) = 2;
            f(1,1) = 10;
            test(!f.isSubmodular());
        }
    }
};

struct GraphicalModelTest
{
    opengm::DiscreteSpace *gm0Space, *gm1Space, *gm2Space; 
    // test space, space for test graphical models
    std::vector<opengm::ExplicitFactor<float> > gm0Potentials, gm1Potentials, gm2Potentials; 
    // potentials of test graphical models
    opengm::GraphicalModel<opengm::ExplicitFactor<float>,opengm::Multiplier > *gm0, *gm1, *gm2; 
    // test graphical models

    void constructionQueryAndAccess()
    {
        // function call
        opengm::GraphicalModel<opengm::ExplicitFactor<float>, opengm::Multiplier> gm1;

        // post condition test
        test(gm1.numberOfFactors() == 0);

        // function call
        size_t nos[] = {2,3,2};
        opengm::DiscreteSpace space1(&nos[0], &nos[3]);
        std::vector<opengm::ExplicitFactor<float> > factors(2);
        size_t vis[] = {0,1};
        factors[0] = opengm::ExplicitFactor<float>(space1, &vis[0], &vis[2]);
        factors[1] = opengm::ExplicitFactor<float>(space1, 2);
        opengm::GraphicalModel<opengm::ExplicitFactor<float>, opengm::Multiplier> gm2(factors.begin(), factors.end());

        // post condition test
        test(&gm2.space() == &space1);
        test(gm2.numberOfFactors() == 2);
        ExplicitFactorTest::testEqualExplicitFactor(gm2[0], factors[0]);
        ExplicitFactorTest::testEqualExplicitFactor(gm2[1], factors[1]);
    }

    void manipulation()
    {
        // function call addFactor
        size_t nos[] = {2,2,2,2};
        opengm::DiscreteSpace space1(&nos[0], &nos[4]);
        opengm::GraphicalModel<opengm::ExplicitFactor<float>, opengm::Multiplier> gm;
        std::vector<opengm::ExplicitFactor<float> > factors(7);
        for(size_t j=0; j<3; ++j) {
            size_t variables[2] = {j, j+1};
            factors[j].assign(space1, &variables[0], &variables[2]);
            size_t k = gm.addFactor(factors[j]);
            test(j == k);
        }
        for(size_t j=0; j<4; ++j) {
            size_t k = 3+j;
            size_t variables[1] = {j};
            factors[k].assign(space1, &variables[0], &variables[1]);
            size_t m = gm.addFactor(factors[k]);
            test(k == m);
        }

        // post condition test
        for(size_t j=0; j<7; ++j) {
            ExplicitFactorTest::testEqualExplicitFactor(gm[j], factors[j]);
        }
    }

    void evaluate()
    {
        //function evaluate()
        size_t nos[] = {2,2,2};
        opengm::DiscreteSpace space1(&nos[0], &nos[3]);
        std::vector<opengm::ExplicitFactor<float> > factors(2);
        size_t vis1[] = {0,1};
        size_t vis2[] = {2};
        factors[0] = opengm::ExplicitFactor<float>(space1, &vis1[0], &vis1[2]);
        factors[1] = opengm::ExplicitFactor<float>(space1, &vis2[0], &vis2[1]);
        factors[0](0,0) = 1;
        factors[0](0,1) = 2;
        factors[0](1,0) = 3;
        factors[0](1,1) = 4;
        factors[1](0)   = 7;
        factors[1](1)   = 9;
        opengm::GraphicalModel<opengm::ExplicitFactor<float>, opengm::Multiplier> gm1(factors.begin(), factors.end());
        opengm::GraphicalModel<opengm::ExplicitFactor<float>, opengm::Adder> gm2(factors.begin(), factors.end());

        std::vector<size_t> x( 3, 0);
        std::vector<size_t> y( 3, 1);
        //float vv = gm1.evaluate(y);
        //test(gm1.evaluate(x) == 1);

        test(gm1.evaluate(x) == 7);
        test(gm1.evaluate(y) == 36);  
        test(gm2.evaluate(x) == 8);
        test(gm2.evaluate(y) == 13);
    }

    void isAcyclic() {
        typedef opengm::DiscreteSpace Space;
        typedef opengm::ExplicitFactor<float> Factor;
        typedef opengm::GraphicalModel<Factor, opengm::Adder> GraphicalModel;
        {
            size_t nos[] = {2, 2, 2, 2};
            Space space(nos, nos+4);
            GraphicalModel gm;
            {
                size_t vi[] = {0,1,2};
                Factor f(space, vi, vi+3);
                gm.addFactor(f);
            }
            {
                size_t vi = 2;
                Factor f(space, &vi, &vi+1);
                gm.addFactor(f);
            }
            {
                size_t vi[] = {2,3};
                Factor f(space, vi, vi+2);
                gm.addFactor(f);
            }
            test(gm.isAcyclic());
        }
        {
            size_t nos[] = {2, 2, 2};
            Space space(nos, nos+3);
            GraphicalModel gm;
            {
                size_t vi[] = {0,1,2};
                Factor f(space, vi, vi+3);
                gm.addFactor(f);
            }
            {
                size_t vi[] = {1,2};
                Factor f(space, vi, vi+2);
                gm.addFactor(f);
            }
            test(!gm.isAcyclic());
        }
        {
            size_t nos[] = {2, 2, 2};
            Space space(nos, nos+3);
            GraphicalModel gm;
            {
                size_t vi[] = {0,1};
                Factor f(space, vi, vi+2);
                gm.addFactor(f);
            }
            {
                size_t vi[] = {1,2};
                Factor f(space, vi, vi+2);
                gm.addFactor(f);
            }
            {
                size_t vi[] = {0,2};
                Factor f(space, vi, vi+2);
                gm.addFactor(f);
            }
            test(!gm.isAcyclic());
            
        }
    }
};


struct BeliefPropagationTest
{
    opengm::DiscreteSpace *gm0Space, *gm1Space, *gm2Space; 
    std::vector<opengm::ExplicitFactor<float> > gm0Potentials, gm1Potentials, gm2Potentials; 
    opengm::GraphicalModel<opengm::ExplicitFactor<float>, opengm::Multiplier> *gm0, *gm1, *gm2; 

    BeliefPropagationTest()
        : gm0Potentials(std::vector<opengm::ExplicitFactor<float> >(3)),
        gm1Potentials(std::vector<opengm::ExplicitFactor<float> >(3)),
        gm2Potentials(std::vector<opengm::ExplicitFactor<float> >(2))
    {
        std::set<size_t> variableIndices;
        std::vector<size_t> numbersOfStates;

        // test graphical models:
        // gm0, tree like:
        numbersOfStates.resize(4);
        numbersOfStates[0] = 2;
        numbersOfStates[1] = 3;
        numbersOfStates[2] = 2;
        numbersOfStates[3] = 2;
        gm0Space = new opengm::DiscreteSpace(numbersOfStates.begin(), numbersOfStates.end());

        variableIndices.insert(0);
        variableIndices.insert(1);
        variableIndices.insert(2);
        gm0Potentials[0] = opengm::ExplicitFactor<float>(*gm0Space, variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        gm0Potentials[0](0,0,0) = 2;
        gm0Potentials[0](0,0,1) = 3;
        gm0Potentials[0](0,1,0) = 4;
        gm0Potentials[0](0,1,1) = 5;
        gm0Potentials[0](0,2,0) = 6;
        gm0Potentials[0](0,2,1) = 7;
        gm0Potentials[0](1,0,0) = 8;
        gm0Potentials[0](1,0,1) = 9;
        gm0Potentials[0](1,1,0) = 10;
        gm0Potentials[0](1,1,1) = 11;
        gm0Potentials[0](1,2,0) = 12;
        gm0Potentials[0](1,2,1) = 13;

        variableIndices.insert(1);
        gm0Potentials[1] = opengm::ExplicitFactor<float>(*gm0Space, variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        gm0Potentials[1](0) = 2;
        gm0Potentials[1](1) = 3;
        gm0Potentials[1](2) = 4;

        variableIndices.insert(2);
        variableIndices.insert(3);
        gm0Potentials[2] = opengm::ExplicitFactor<float>(*gm0Space, variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        gm0Potentials[2](0,0) = 2;
        gm0Potentials[2](0,1) = 3;
        gm0Potentials[2](1,0) = 4;
        gm0Potentials[2](1,1) = 5;

        gm0 = new opengm::GraphicalModel<opengm::ExplicitFactor<float>, opengm::Multiplier>(gm0Potentials.begin(),gm0Potentials.end());

        // graphical model 1, loopy:
        numbersOfStates.resize(3);
        numbersOfStates[0] = 2;
        numbersOfStates[1] = 3;
        numbersOfStates[2] = 2;
        gm1Space = new opengm::DiscreteSpace(numbersOfStates.begin(), numbersOfStates.end());

        variableIndices.insert(0);
        variableIndices.insert(1);
        variableIndices.insert(2);
        gm1Potentials[0] = opengm::ExplicitFactor<float>(*gm1Space, variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        gm1Potentials[0](0,0,0) = 2;
        gm1Potentials[0](0,0,1) = 3;
        gm1Potentials[0](0,1,0) = 4;
        gm1Potentials[0](0,1,1) = 5;
        gm1Potentials[0](0,2,0) = 6;
        gm1Potentials[0](0,2,1) = 7;
        gm1Potentials[0](1,0,0) = 8;
        gm1Potentials[0](1,0,1) = 9;
        gm1Potentials[0](1,1,0) = 10;
        gm1Potentials[0](1,1,1) = 11;
        gm1Potentials[0](1,2,0) = 12;
        gm1Potentials[0](1,2,1) = 13;
        gm1Potentials[0].normalize();

        variableIndices.insert(0);
        gm1Potentials[1] = opengm::ExplicitFactor<float>(*gm1Space, variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        gm1Potentials[1](0) = 2;
        gm1Potentials[1](1) = 3;

        variableIndices.insert(1);
        variableIndices.insert(2);
        gm1Potentials[2] = opengm::ExplicitFactor<float>(*gm1Space, variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        gm1Potentials[2](0,0) = 2;
        gm1Potentials[2](0,1) = 3;
        gm1Potentials[2](1,0) = 4;
        gm1Potentials[2](1,1) = 5;
        gm1Potentials[2](2,0) = 6;
        gm1Potentials[2](2,1) = 7;
        gm1Potentials[2].normalize();

        gm1 = new opengm::GraphicalModel<opengm::ExplicitFactor<float>, opengm::Multiplier>(gm1Potentials.begin(),gm1Potentials.end());

        // graphical model 2, degenerate:
        numbersOfStates.resize(2);
        numbersOfStates[0] = 1;
        numbersOfStates[1] = 1;
        gm2Space = new opengm::DiscreteSpace(numbersOfStates.begin(), numbersOfStates.end());

        gm2Potentials[0] = opengm::ExplicitFactor<float>(*gm2Space, 2); // constant potential

        variableIndices.insert(1);
        gm2Potentials[1] = opengm::ExplicitFactor<float>(*gm2Space, variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        gm2Potentials[1](0) = 2;

        gm2 = new opengm::GraphicalModel<opengm::ExplicitFactor<float>, opengm::Multiplier>(gm2Potentials.begin(),gm2Potentials.end());
    }

    ~BeliefPropagationTest() {
        delete gm0;
        delete gm1;
        delete gm2;

        delete gm0Space;
        delete gm1Space;
        delete gm2Space;
    }

    void initializeGm0(std::vector<opengm::ExplicitFactor<float> >& mxf, std::vector<opengm::ExplicitFactor<float> >& mfx)
    {
        mxf.resize(6);
        for(size_t j=0; j<6; ++j) {
            mxf[j].assign(gm0->space(), 1);
        }   

        mfx.resize(6);
        std::set<size_t> variableIndices;
        variableIndices.insert(0);
        mfx[0].assign(gm0->space(), variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        variableIndices.insert(1);
        mfx[1].assign(gm0->space(), variableIndices.begin(), variableIndices.end());
        mfx[3].assign(gm0->space(), variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        variableIndices.insert(2);
        mfx[2].assign(gm0->space(), variableIndices.begin(), variableIndices.end());
        mfx[4].assign(gm0->space(), variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        variableIndices.insert(3);
        mfx[5].assign(gm0->space(), variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        for(size_t j=0; j<6; ++j) {
            mfx[j].normalize();
        }   
    }

    // auxiliary functions for belief propagation test

    void propagateGm0
        (
        std::vector<opengm::ExplicitFactor<float> >& mxf,
        std::vector<opengm::ExplicitFactor<float> >& mfx
        ) 
    {
        std::vector<opengm::ExplicitFactor<float> > mxft(6), mfxt(6);
        std::set<size_t> indices;

        // mxf update
        mxft[0].assign(gm0->space(), 1);
        mxft[1] = mfx[3];
        mxft[2] = mfx[4];
        mxft[3] = mfx[1];
        mxft[4] = mfx[2];
        mxft[5].assign(gm0->space(), 1);

        // mfx update:
        mfxt[0] = gm0Potentials[0] * mxf[1] * mxf[2];
        indices.insert(1);
        indices.insert(2);  
        mfxt[0].maximize(indices.begin(), indices.end());
        indices.clear();
        mfxt[0].normalize();

        mfxt[1] = gm0Potentials[0] * mxf[0] * mxf[2];
        indices.insert(0);
        indices.insert(2);
        mfxt[1].maximize(indices.begin(), indices.end());
        indices.clear();
        mfxt[1].normalize();

        mfxt[2] = gm0Potentials[0] * mxf[0] * mxf[1];
        indices.insert(0);
        indices.insert(1);
        mfxt[2].maximize(indices.begin(), indices.end());
        indices.clear();
        mfxt[2].normalize();

        mfxt[3] = gm0Potentials[1];
        mfxt[3].normalize();

        mfxt[4] = gm0Potentials[2] * mxf[5];
        indices.insert(3);
        mfxt[4].maximize(indices.begin(), indices.end());
        indices.clear();
        mfxt[4].normalize();

        mfxt[5] = gm0Potentials[2] * mxf[4];
        indices.insert(2);
        mfxt[5].maximize(indices.begin(), indices.end());
        indices.clear();
        mfxt[5].normalize();

        // output
        for(size_t j=0; j<6; ++j) {
            mxf[j] = mxft[j];
            mfx[j] = mfxt[j];
        }
    }

    void beliefGm0(std::vector<opengm::ExplicitFactor<float> >& belief, std::vector<opengm::ExplicitFactor<float> >& mfx)
    {
        belief[0] = mfx[0];
        belief[1] = (mfx[1] * mfx[3]);
        belief[1].normalize();
        belief[2] = (mfx[2] * mfx[4]);
        belief[2].normalize();
        belief[3] = mfx[5];
    }

    void initializeGm1(std::vector<opengm::ExplicitFactor<float> >& mxf, std::vector<opengm::ExplicitFactor<float> >& mfx)
    {
        mxf.resize(6);
        for(size_t j=0; j<6; ++j){
            mxf[j].assign(gm1->space(), 1);
        }   

        mfx.resize(6);
        std::set<size_t> variableIndices;
        variableIndices.insert(0);
        mfx[0].assign(gm1->space(), variableIndices.begin(), variableIndices.end());
        mfx[3].assign(gm1->space(), variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        variableIndices.insert(1);
        mfx[1].assign(gm1->space(), variableIndices.begin(), variableIndices.end());
        mfx[4].assign(gm1->space(), variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        variableIndices.insert(2);
        mfx[2].assign(gm1->space(), variableIndices.begin(), variableIndices.end());
        mfx[5].assign(gm1->space(), variableIndices.begin(), variableIndices.end());
        variableIndices.clear();
        for(size_t j=0; j<6; ++j) {
            mfx[j].normalize();
        }   
    }

    void propagateGm1(std::vector<opengm::ExplicitFactor<float> >& mxf, std::vector<opengm::ExplicitFactor<float> >& mfx)
    {
        std::vector<opengm::ExplicitFactor<float> > mxft(6), mfxt(6);
        std::set<size_t> indices;

        // mxf update
        mxft[0] = mfx[3];
        mxft[1] = mfx[4];
        mxft[2] = mfx[5];
        mxft[3] = mfx[0];
        mxft[4] = mfx[1];
        mxft[5] = mfx[2];

        // mfx update:
        mfxt[0] = gm1Potentials[0] * mxf[1] * mxf[2];
        indices.insert(1);
        indices.insert(2);
        mfxt[0].maximize(indices.begin(), indices.end());
        indices.clear();
        mfxt[0].normalize();

        mfxt[1] = gm1Potentials[0] * mxf[0] * mxf[2];
        indices.insert(0);
        indices.insert(2);
        mfxt[1].maximize(indices.begin(), indices.end());
        indices.clear();
        mfxt[1].normalize();

        mfxt[2] = gm1Potentials[0] * mxf[0] * mxf[1];
        indices.insert(0);
        indices.insert(1);
        mfxt[2].maximize(indices.begin(), indices.end());
        indices.clear();
        mfxt[2].normalize();

        mfxt[3] = gm1Potentials[1];
        mfxt[3].normalize();

        mfxt[4] = gm1Potentials[2] * mxf[5];
        indices.insert(2);
        mfxt[4].maximize(indices.begin(), indices.end());
        indices.clear();
        mfxt[4].normalize();

        mfxt[5] = gm1Potentials[2] * mxf[4];
        indices.insert(1);
        mfxt[5].maximize(indices.begin(), indices.end());
        indices.clear();
        mfxt[5].normalize();

        // output
        for(size_t j=0; j<6; ++j) {
            mxf[j] = mxft[j];
            mfx[j] = mfxt[j];
        }
    }

    void beliefGm1(std::vector<opengm::ExplicitFactor<float> >& belief, std::vector<opengm::ExplicitFactor<float> >& mfx)
    {
        belief[0] = (mfx[0] * mfx[3]);
        belief[0].normalize();
        belief[1] = (mfx[1] * mfx[4]);
        belief[1].normalize();
        belief[2] = (mfx[2] * mfx[5]);
        belief[2].normalize();
    }

    void initializeGm2(std::vector<opengm::ExplicitFactor<float> >& mxf, std::vector<opengm::ExplicitFactor<float> >& mfx)
    {
        mxf.resize(1);
        mxf[0].assign(gm2->space(), 1);

        mfx.resize(1);
        std::set<size_t> variableIndices;
        variableIndices.insert(1);
        mfx[0].assign(gm2->space(), variableIndices.begin(), variableIndices.end());
        mfx[0].normalize();
    }

    void propagateGm2(std::vector<opengm::ExplicitFactor<float> >& mxf, std::vector<opengm::ExplicitFactor<float> >& mfx)
    {
        mxf[0].assign(gm2->space(), 1);

        mfx[0] = gm2Potentials[1];
        mfx[0].normalize();
    }

    void beliefGm2(std::vector<opengm::ExplicitFactor<float> >& belief, std::vector<opengm::ExplicitFactor<float> >& mfx)
    {
        std::set<size_t> variableIndices;
        variableIndices.insert(0);
        belief[0] = opengm::ExplicitFactor<float>(gm2->space(), variableIndices.begin(), variableIndices.end());
        belief[0].normalize();
        belief[1] = mfx[0];
    }

    // end of: auxiliary functions for belief propagation test

    void beliefPropagationTest(void)
    {
        opengm::ExplicitFactor<float> belief; // computed belief
        std::vector<opengm::ExplicitFactor<float> > trueBelief; // true belief
        std::vector<opengm::ExplicitFactor<float> > mxf, mfx; // true messages
        std::vector<size_t> vi;
        // gm0:
        trueBelief.resize(4);
        initializeGm0(mxf, mfx);
        opengm::BeliefPropagation<opengm::GraphicalModel<opengm::ExplicitFactor<float>,opengm::Multiplier>,
            opengm::Maximizer,
            opengm::MaxDistance> 
            bp0(*gm0); 

        for(size_t step=0; step<25; ++step) {
            // compare:
            beliefGm0(trueBelief, mfx);
            for(size_t var=0; var<gm0->space().dimension(); ++var) {
                bp0.marginal(var, belief);
                test(trueBelief[var].numberOfVariables() == belief.numberOfVariables());
                vi.resize(trueBelief[var].numberOfVariables());
                trueBelief[var].variableIndices(vi.begin());
                std::vector<size_t> vi2(belief.numberOfVariables());
                belief.variableIndices(vi2.begin());
                testEqualSequence(vi.begin(), vi.end(), vi2.begin());
                for(size_t j=0; j<gm0->space().numberOfStates(var); ++j){
                    testEqualTolerance(trueBelief[var](j), belief(j), 1e-6);
                }
            }
            // propagate:
            propagateGm0(mxf, mfx); // truth 
            bp0.propagate();
            float dump = bp0.convergence();
            DGM_UNUSED(dump);
        } 

        // gm1:
        trueBelief.resize(3);
        initializeGm1(mxf, mfx);
        opengm::BeliefPropagation<opengm::GraphicalModel<opengm::ExplicitFactor<float>,opengm::Multiplier>,
            opengm::Maximizer,
            opengm::MaxDistance> 
            bp1(*gm1);
        for(size_t step=0; step<25; ++step) {
            // compare:
            beliefGm1(trueBelief, mfx);
            for(size_t var=0; var<gm1->space().dimension(); ++var) {
                bp1.marginal(var, belief);
                test(trueBelief[var].numberOfVariables() == belief.numberOfVariables());
                vi.resize(trueBelief[var].numberOfVariables());
                trueBelief[var].variableIndices(vi.begin());
                std::vector<size_t> vi2(belief.numberOfVariables());
                belief.variableIndices(vi2.begin());
                testEqualSequence(vi.begin(), vi.end(), vi2.begin());
                for(size_t j=0; j<gm1->space().numberOfStates(var); ++j) {
                    testEqualTolerance(trueBelief[var](j), belief(j), 1e-6);
                }
            }
            // propagate:
            propagateGm1(mxf, mfx); // truth
            bp1.propagate();
            float dump = bp1.convergence();
            DGM_UNUSED(dump);
        }

        // gm2:
        trueBelief.resize(2);
        initializeGm2(mxf, mfx);
        opengm::BeliefPropagation<opengm::GraphicalModel<opengm::ExplicitFactor<float>,opengm::Multiplier>,
            opengm::Maximizer,
            opengm::MaxDistance> 
            bp2(*gm2);
        for(size_t step=0; step<25; ++step) {
            // compare:
            beliefGm2(trueBelief, mfx);
            for(size_t var=0; var<gm2->space().dimension(); ++var) {
                bp2.marginal(var, belief);
                test(trueBelief[var].numberOfVariables() == belief.numberOfVariables());
                vi.resize(trueBelief[var].numberOfVariables());
                trueBelief[var].variableIndices(vi.begin());
                std::vector<size_t> vi2(belief.numberOfVariables());
                belief.variableIndices(vi2.begin());
                testEqualSequence(vi.begin(), vi.end(), vi2.begin());
                for(size_t j=0; j<gm2->space().numberOfStates(var); ++j) {
                    testEqual(trueBelief[var](j), belief(j));
                }
            }
            // propagate:
            propagateGm2(mxf, mfx); // truth
            bp2.propagate();
            float dump = bp2.convergence();
            DGM_UNUSED(dump);
        }
    }
};

struct TreeReweightedBeliefPropagationTest
{
    void run()
    {
        size_t nos[] = {2,2,2};
	opengm::DiscreteSpace space1(&nos[0], &nos[3]);
	std::vector<opengm::ExplicitFactor<float> > factors(5);
	size_t vis1[] = {0,1};
	size_t vis2[] = {1,2};
	size_t vis3[] = {0,2};
	size_t vis4[] = {0};
	size_t vis5[] = {0,1,2};
	factors[0] = opengm::ExplicitFactor<float>(space1, &vis1[0], &vis1[2]);
	factors[1] = opengm::ExplicitFactor<float>(space1, &vis2[0], &vis2[2]);
       	factors[2] = opengm::ExplicitFactor<float>(space1, &vis3[0], &vis3[2]);
      	factors[3] = opengm::ExplicitFactor<float>(space1, &vis4[0], &vis4[1]);
      	factors[4] = opengm::ExplicitFactor<float>(space1, &vis5[0], &vis5[3]);
        factors[0](0,0) = 1;
        factors[0](0,1) = 0;
        factors[0](1,0) = 0;
        factors[0](1,1) = 1;
        factors[1](0,0) = 1;
        factors[1](0,1) = 0;
        factors[1](1,0) = 0;
        factors[1](1,1) = 1;
        factors[2](0,0) = 1;
        factors[2](0,1) = 0;
        factors[2](1,0) = 0;
        factors[2](1,1) = 1;
        factors[3](0) = 0;
        factors[3](1) = 0;
	factors[4](0,0,0) = 1;
        factors[4](0,0,1) = 0;
        factors[4](0,1,0) = 1;
        factors[4](0,1,1) = 1;
	factors[4](1,0,0) = 1;
        factors[4](1,0,1) = 1;
        factors[4](1,1,0) = 1;
        factors[4](1,1,1) = 1;
        opengm::GraphicalModel<opengm::ExplicitFactor<float>,opengm::Multiplier> gm1(factors.begin(), factors.end());
        opengm::GraphicalModel<opengm::ExplicitFactor<float>,opengm::Adder>      gm2(factors.begin(), factors.end());

	std::vector<size_t> sol;
	typedef opengm::TreeReweightedBeliefPropagation<
	  opengm::GraphicalModel<opengm::ExplicitFactor<float>,opengm::Adder>, 
	  opengm::Minimizer,
	  opengm::MaxDistance
	> TRBP;
	typedef TRBP::Parameter Para;

	Para para;
	para.rho_.resize(5);
	para.rho_[0]=2.0/4;
	para.rho_[1]=2.0/4;
	para.rho_[2]=2.0/4;
	para.rho_[3]=1;
	para.rho_[4]=1.0/4;
	TRBP trbp(gm2,para);
	trbp.infer();
	test(trbp.arg(sol)==opengm::NORMAL);
	//test(sol[0]==0 && sol[1]==0 && sol[2]==1 );
	float bound;
	trbp.bound(bound);
	test(bound<=1); 
    }
};

struct AStarTest
{
    void run()
    {
        //function evaluate()
        size_t nos[] = {2,2,2};
        opengm::DiscreteSpace space1(&nos[0], &nos[3]);
        std::vector<opengm::ExplicitFactor<float> > factors(2);
        size_t vis1[] = {0,1};
        size_t vis2[] = {2};
        factors[0] = opengm::ExplicitFactor<float>(space1, &vis1[0], &vis1[2]);
        factors[1] = opengm::ExplicitFactor<float>(space1, &vis2[0], &vis2[1]);
        factors[0](0,0) = 1;
        factors[0](0,1) = 2;
        factors[0](1,0) = 3;
        factors[0](1,1) = 4;
        factors[1](0)   = 7;
        factors[1](1)   = 9;
        opengm::GraphicalModel<opengm::ExplicitFactor<float>,opengm::Multiplier> gm1(factors.begin(), factors.end());
        opengm::GraphicalModel<opengm::ExplicitFactor<float>,opengm::Adder>      gm2(factors.begin(), factors.end());

        std::vector<size_t> sol;    
        opengm::AStar<opengm::GraphicalModel<opengm::ExplicitFactor<float>,opengm::Multiplier>, opengm::Maximizer> astar1(gm1); 
        astar1.infer();
        test(astar1.arg(sol)==opengm::NORMAL);
        test(sol.size()==3);
        test(sol[0]==1 && sol[1]==1 && sol[2]==1 ); 

        opengm::AStar<opengm::GraphicalModel<opengm::ExplicitFactor<float>,opengm::Multiplier>, opengm::Minimizer> astar2(gm1); 
        astar2.infer();
        test(astar2.arg(sol)==opengm::NORMAL);   
        test(sol.size()==3);
        test(sol[0]==0 && sol[1]==0 && sol[2]==0 );
    }
};

struct MovemakerTest {
    void run() {
        // space
        size_t numbersOfStates[] = {3, 3, 3, 3, 3};
        opengm::DiscreteSpace space(numbersOfStates, numbersOfStates+5);

        // graphical model
        typedef opengm::ExplicitFactor<float> Factor;
        typedef opengm::GraphicalModel<Factor, opengm::Adder> GraphicalModel;
        GraphicalModel gm;

        // single site factors
        for(size_t j=0; j<space.dimension(); ++j) {
            size_t variableIndices[] = {j};
            Factor f(space, variableIndices, variableIndices+1);
            f(0) = 0.0f;
            f(1) = 0.2f;
            f(2) = 0.2f;
            gm.addFactor(f);
        }

        // 2nd order factors
        for(size_t j=0; j<space.dimension()-1; ++j) {
            size_t variableIndices[] = {j, j+1};
            Factor f(space, variableIndices, variableIndices+2);
            for(size_t j=0; j<9; ++j) {
                f(j) = static_cast<float>(j) / 10;
            }
            gm.addFactor(f);
        }

        // 3rd order factor
        {
            size_t variableIndices[] = {0, 2, 4};
            Factor f(space, variableIndices, variableIndices+3, 1.0f);
            for(size_t j=0; j<27; ++j) {
                f(j) = static_cast<float>(j) / 20;
            }
            gm.addFactor(f);
        }

        // Move maker test by exhaustive comparison
        typedef opengm::Movemaker<GraphicalModel> Movemaker;
        Movemaker movemaker(gm);
        std::vector<size_t> state(gm.space().dimension());
        std::vector<size_t> vi(gm.space().dimension());
        for(size_t j=0; j<gm.space().dimension(); ++j) {
            vi[j] = j;
        }
        bool overflow = false;
        while(!overflow) {
            test(
                movemaker.energyAfterMove(vi.begin(), vi.end(), state.begin())
                == gm.evaluate(state.begin()) 
            );
            for(size_t j=0; j<gm.space().dimension(); ++j) {
                if(state[j]+1 < gm.space().numberOfStates(j)) {
                    ++state[j];
                    break;
                }
                else {
                    state[j] = 0;
                    if(j == gm.space().dimension()-1) {
                        overflow = true;
                    }
                }
            }
        }
    }
};

class LazyFlipperTest {
public:
    typedef double Energy;
    typedef opengm::DiscreteSpace Space;
    typedef opengm::ExplicitFactor<Energy> Factor;
    typedef opengm::GraphicalModel<Factor, opengm::Adder> GraphicalModel;
    typedef opengm::LazyFlipper<GraphicalModel, opengm::Minimizer> LazyFlipper;

    LazyFlipperTest() {
        // build the discrete space
        std::vector<size_t> numbersOfStates(6, 2);
        space_ = Space(numbersOfStates.begin(), numbersOfStates.end());

        // single site factors
        for(size_t j=0; j<space_.dimension(); ++j) {
            Factor factor(space_, &j, &j + 1);
            factor(0) = static_cast<Energy>(rand()) / static_cast<Energy>(RAND_MAX);
            factor(1) = 1.0 - factor(0);
            gm_.addFactor(factor);
        }

        // Potts factors
        Energy alpha = 0.2f;
        for(size_t j=0; j<2; ++j) {
            {
                size_t vi[] = {j, j+1};
                Factor factor(space_, vi, vi+2);
                factor(0,0) = 0.0f;
                factor(0,1) = alpha;
                factor(1,0) = alpha;
                factor(1,1) = 0.0f;
                gm_.addFactor(factor);
            }
            {
                size_t vi[] = {j+3, j+4};
                Factor factor(space_, vi, vi+2);
                factor(0,0) = 0.0f;
                factor(0,1) = alpha;
                factor(1,0) = alpha;
                factor(1,1) = 0.0f;
                gm_.addFactor(factor);
            }
        }
        for(size_t j=0; j<3; ++j) {
            size_t vi[] = {j, j+3};
            Factor factor(space_, vi, vi+2);
            factor(0,0) = 0.0f;
            factor(0,1) = alpha;
            factor(1,0) = alpha;
            factor(1,1) = 0.0f;
            gm_.addFactor(factor);
        }
    }

    void run() {
        LazyFlipper flipper(gm_);
        size_t maxSubgraphSize = 6;
        flipper.setMaxSubgraphSize(maxSubgraphSize);
        opengm::LazyFlipperVisitor<LazyFlipper> visitor;
        flipper.infer(visitor);

        std::vector<GraphicalModel::state_type> state(gm_.space().dimension());
        flipper.arg(state);

        Energy energy = 0;
        energy = flipper.optimum();
    }

private:
    GraphicalModel gm_;
    Space space_;
};

struct ICMTest
{
    typedef float value_type;
    typedef opengm::ExplicitFactor<value_type> Factor;
    typedef opengm::DiscreteSpace Space;

    void run()
    {
        //function evaluate()
        size_t nos[] = {2,2,2};
        Space space1(&nos[0], &nos[3]);
        std::vector<Factor> factors(2);
        size_t vis1[] = {0,1};
        size_t vis2[] = {2};
        factors[0] = Factor(space1, &vis1[0], &vis1[2]);
        factors[1] = Factor(space1, &vis2[0], &vis2[1]);
        factors[0](0,0) = 1;
        factors[0](0,1) = 2;
        factors[0](1,0) = 3;
        factors[0](1,1) = 4;
        factors[1](0)   = 7;
        factors[1](1)   = 9;
        {
            typedef opengm::GraphicalModel<Factor, opengm::Multiplier> GraphicalModel;
            typedef opengm::ICM<GraphicalModel, opengm::Maximizer> ICM;
            GraphicalModel gm1(factors.begin(), factors.end());
            ICM icm1(gm1); 
            icm1.infer();
            std::vector<size_t> sol;
            test(icm1.arg(sol) == opengm::NORMAL); 
            test(sol[0]==1 && sol[1]==1 && sol[2]==1);
        }
    }
};

template<class FACTOR>
class SerializationTest {
public:
    typedef FACTOR Factor;
    typedef typename Factor::value_type value_type;
    typedef opengm::DiscreteSpace Space;
    typedef opengm::GraphicalModel<Factor, opengm::Adder> GraphicalModel;

    SerializationTest()
    : gm_(GraphicalModel(4)) {
        size_t numbersOfStates[] = {2, 2, 2, 2};
        space_ = Space(numbersOfStates, numbersOfStates+4);

        // single site factors
        {
            marray::Vector<value_type> table(2);
            table(0) = 1;
            table(1) = 2;
            for(size_t j=0; j<3; ++j) {
                Factor f(space_, &j, &j+1, table.begin(), table.end());
                gm_.addFactor(f);
            }
        }

        // 2nd order factors
        {
            marray::Matrix<value_type> table(2, 2);
            table(0,0) = 1;
            table(1,0) = 2;
            table(0,1) = 3;
            table(1,1) = 4;
            size_t vi[2];
            for(size_t j=0; j<3; ++j) {
                vi[0] = j;
                vi[1] = j+1;
                Factor f(space_, vi, vi+2, table.begin(), table.end());
                gm_.addFactor(f);
            }
        }
    }

    void run() {
        // serialize
        marray::Vector<value_type> serialization;
        opengm::serialize(gm_, serialization);

        // de-serialize
        Space space;
        GraphicalModel gm;
        opengm::deserialize(serialization, gm, space);

        // compare
        test(gm.space().dimension() == gm_.space().dimension());
        test(gm.numberOfFactors() == gm_.numberOfFactors());
        for(size_t j=0; j<gm.numberOfFactors(); ++j) {
            const Factor& f1 = gm[j];
            const Factor& f2 = gm_[j];
            test(f1.numberOfVariables() == f2.numberOfVariables());
            for(size_t k=0; k<f1.numberOfVariables(); ++k) {
                test(f1.variableIndex(k) == f2.variableIndex(k));
            }
            test(f1.table().size() == f2.table().size());
            for(size_t k=0; k<f1.table().size(); ++k) {
                test(f1(k) == f2(k));
            }
        }
    }

private:
    Space space_;
    GraphicalModel gm_;
};

class DecomposerTest {
public:
    typedef double value_type;
    typedef opengm::ExplicitFactor<value_type> Factor;
    typedef opengm::DiscreteSpace Space;
    typedef opengm::GraphicalModel<Factor, opengm::Adder> GraphicalModel;

    DecomposerTest()
    {
        size_t numbersOfStates[] = {2,2,2,2,2,2};
        space_ = Space(numbersOfStates, numbersOfStates+6);

        // unary factors
        for(size_t j=0; j<space_.dimension(); ++j) {
            value_type values[] = {0.4, 0.6};
            Factor f(space_, &j, &j+1, values, values+2);
            gm_.addFactor(f);
        }

        // 2nd order factors
        value_type values[] = {0.2, 0.3, 0.4, 0.1};
        size_t vi[2];
        { 
            vi[0] = 0; vi[1] = 1; 
            Factor f(space_, vi, vi+2, values, values+4);
            gm_.addFactor(f); 
        }
        { 
            vi[0] = 0; vi[1] = 3; 
            Factor f(space_, vi, vi+2, values, values+4);
            gm_.addFactor(f); 
        }
        { 
            vi[0] = 1; vi[1] = 2; 
            Factor f(space_, vi, vi+2, values, values+4);
            gm_.addFactor(f); 
        }
        { 
            vi[0] = 1; vi[1] = 4; 
            Factor f(space_, vi, vi+2, values, values+4);
            gm_.addFactor(f); 
        }
        { 
            vi[0] = 2; vi[1] = 5; 
            Factor f(space_, vi, vi+2, values, values+4);
            gm_.addFactor(f); 
        }
        { 
            vi[0] = 3; vi[1] = 4; 
            Factor f(space_, vi, vi+2, values, values+4);
            gm_.addFactor(f); 
        }
        { 
            vi[0] = 4; vi[1] = 5; 
            Factor f(space_, vi, vi+2, values, values+4);
            gm_.addFactor(f); 
        }
    }

    void run() {
        std::vector<std::vector<size_t> > decomposition;
        opengm::Decomposer<GraphicalModel> decomposer;
        decomposer.decompose(gm_, decomposition);
        
        /*
        std::cout << "decomposition into " << decomposition.size() << " acyclic subgraphs." << std::endl;
        for(size_t j=0; j<decomposition.size(); ++j) {
            std::cout << "subgraph " << j << std::endl;
            for(size_t k=0; k<decomposition[j].size(); ++k) {
                std::cout << decomposition[j][k] << ' ';
            }
            std::cout << std::endl;
        }
        */

        std::vector<bool> foundFactor(gm_.numberOfFactors());
        for(size_t j=0; j<decomposition.size(); ++j) {
            for(size_t k=0; k<decomposition[j].size(); ++k) {
                foundFactor[decomposition[j][k]] = true;
            }
        }
        for(size_t j=0; j<foundFactor.size(); ++j) {
            test(foundFactor[j] == true);
        }
    }

private:
    Space space_;
    GraphicalModel gm_;
};

int main()
{
    std::cout << "ExplicitFactorTest... ";
    { ExplicitFactorTest t; t.constructionAndQuery(); }
    { ExplicitFactorTest t; t.assignment(); }
    { ExplicitFactorTest t; t.manipulation(); }
    { ExplicitFactorTest t; t.operation(); }
    { ExplicitFactorTest t; t.accumulation(); }
    { ExplicitFactorTest t; t.maxMinInt(); }
    { ExplicitFactorTest t; t.normalizeSubtractOffset(); }
    { ExplicitFactorTest t; t.isSubmodularTest(); }
    std::cout << "passed." << std::endl;

    std::cout << "GraphicalModelTest... ";
    { GraphicalModelTest t; t.constructionQueryAndAccess(); }
    { GraphicalModelTest t; t.manipulation(); }
    { GraphicalModelTest t; t.evaluate(); }   
    { GraphicalModelTest t; t.isAcyclic(); }   
    std::cout << "passed." << std::endl;

    std::cout << "SerializationTest... ";
    { SerializationTest<opengm::ExplicitFactor<double> > t; t.run(); }
    //{ SerializationTest<opengm::ExplicitFactor<double, ArrayPow2<double> > > t; t.run(); }
    std::cout << "passed." << std::endl;

    std::cout << "ICMTest... ";
    { ICMTest t; t.run(); }
    std::cout << "passed." << std::endl;

    std::cout << "BeliefPropagationTest... ";
    { BeliefPropagationTest t; t.beliefPropagationTest(); }
    std::cout << "passed." << std::endl;

    std::cout << "TreeReweightedBeliefPropagationTest... ";
    { TreeReweightedBeliefPropagationTest t; t.run(); }
    std::cout << "passed." << std::endl;

    std::cout << "AStarTest... ";
    { AStarTest t; t.run(); } 
    std::cout << "passed." << std::endl;

    std::cout << "MoveMakerTest... ";
    { MovemakerTest t; t.run(); } 
    std::cout << "passed." << std::endl;

    std::cout << "LazyFlipperTest... ";
    { LazyFlipperTest t; t.run(); }
    std::cout << "passed." << std::endl;

    std::cout << "DecomposerTest... ";
    { DecomposerTest t; t.run(); }
    std::cout << "passed." << std::endl;

/*
    size_t shape[] = {2, 2};
    ArrayPow2<double> a(shape, shape+2);
    a(0,0) = 1;
    a(0,1) = 2;
    a(1,0) = 3;
    a(1,1) = 4;
    std::cout << a(0,0) << a(0,1) << a(1,0) << a(1,1) << std::endl;

    double q = 0;
    std::vector<size_t> pos(2);

    pos[0] = 0; pos[1] = 0;
    q = a(pos.begin());
    std::cout << q;
  
    pos[0] = 0; pos[1] = 1;
    q = a(pos.begin());
    std::cout << q;

    pos[0] = 1; pos[1] = 0;
    q = a(pos.begin());
    std::cout << q;

    pos[0] = 1; pos[1] = 1;
    q = a(pos.begin());
    std::cout << q;
*/

    return 0;
}
