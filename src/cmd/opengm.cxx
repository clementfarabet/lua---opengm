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
#include <iostream>
#include <string>

#include <marray/marray_hdf5.hxx>
#include <arraypow2/arraypow2.hxx>

#include <tclap/CmdLine.h>

#include <opengm/graphicalmodel.hxx>
#include <opengm/serialization.hxx>
#include <opengm/adder.hxx>
#include <opengm/minimizer.hxx>
#include <opengm/maxdistance.hxx>
#include <opengm/decomposer.hxx>
#include <opengm/inference/beliefpropagation.hxx>
#include <opengm/inference/treereweightedbeliefpropagation.hxx>
#include <opengm/inference/icm.hxx>
#include <opengm/inference/lazyflipper.hxx>
#include <opengm/inference/astar.hxx>

#ifdef WITH_BOOST
    #include <opengm/inference/graphcut.hxx>
#endif
#include <opengm/timing.hxx>
USETICTOC;

template<class ICM>
class ICMProtocolVisitor {
public:
    typedef ICM icm_type;
    typedef typename icm_type::value_type value_type;

    ICMProtocolVisitor(bool protocolStates = false)
        : protocolStates_(protocolStates)
        { TIC; }
    template<class StateIterator>
    void operator()(icm_type& icm, StateIterator stateBegin, StateIterator stateEnd, const value_type& value)
        {
            double t = TOCN;
            if(runtimes_.size() == 0) {
                runtimes_.push_back(t);
            }
            else {
                runtimes_.push_back(t + runtimes_[runtimes_.size()-1]);
            }
            values_.push_back(value);
            if(protocolStates_) {
                size_t numberOfVariables = stateEnd - stateBegin;
                states_.push_back(std::vector<size_t>(numberOfVariables));
                StateIterator it = stateBegin;
                for(size_t j=0; j<numberOfVariables; ++j, ++it) {
                    states_[states_.size()-1][j] = *it;
                }
            }
            std::cout << "step " << runtimes_.size() 
                << ": E=" <<  value
                << " after " << runtimes_[runtimes_.size()-1]/1e3 << " s." << std::endl;
            TIC;
        }
    const std::vector<double>& runtimes() const
        { return runtimes_; }
    const std::vector<value_type>& values() const
        { return values_; }
    const std::vector<std::vector<size_t> >& states() const
        { return states_; }

private:
    std::vector<double> runtimes_;
    std::vector<value_type> values_;
    bool protocolStates_;
    std::vector<std::vector<size_t> > states_;
};

template<class LAZYFLIPPER>
class LazyFlipperProtocolVisitor {
public:
    typedef LAZYFLIPPER lazy_flipper_type;
    typedef typename lazy_flipper_type::value_type value_type;

    LazyFlipperProtocolVisitor(bool protocolStates = false)
        : protocolStates_(protocolStates)
        { TIC; }
    template<class StateIterator>
    void operator()(const lazy_flipper_type& lazyFlipper,
        StateIterator stateBegin, StateIterator stateEnd,
        const value_type& value, const size_t& subgraphSize,
        const size_t& subgraphForestSize)
        {
            double t = TOCN;
            if(runtimes_.size() == 0) { 
                runtimes_.push_back(t);
            }
            else {
                runtimes_.push_back(t + runtimes_[runtimes_.size()-1]);
            }
            values_.push_back(value);
            subgraphSizes_.push_back(subgraphSize);
            subgraphForestSizes_.push_back(subgraphForestSize);
            if(protocolStates_) {
                size_t numberOfVariables = stateEnd - stateBegin;
                if(subgraphSizes_.size() == 1 || subgraphSizes_[subgraphSizes_.size()-2] != subgraphSize) {
                    states_.push_back(std::vector<size_t>(numberOfVariables));
                }
                StateIterator it = stateBegin;
                for(size_t j=0; j<numberOfVariables; ++j, ++it) {
                    states_[states_.size()-1][j] = *it;
                }
            }
            std::cout << "step " << runtimes_.size() 
                << ": subgraph_size=" << subgraphSize
                << ", cs_tree_size=" << subgraphForestSize
                << ", E=" <<  value
                << " after " << runtimes_[runtimes_.size()-1]/1e3 << " s." << std::endl;
            TIC;
        }
    const std::vector<double>& runtimes() const
        { return runtimes_; }
    const std::vector<value_type>& values() const
        { return values_; }
    const std::vector<size_t>& subgraphSizes() const
        { return subgraphSizes_; }
    const std::vector<size_t>& subgraphForestSizes() const
        { return subgraphForestSizes_; }
    const std::vector<std::vector<size_t> >& states() const
        { return states_; }

private:
    std::vector<double> runtimes_;
    std::vector<value_type> values_;
    std::vector<size_t> subgraphSizes_;
    std::vector<size_t> subgraphForestSizes_;
    bool protocolStates_;
    std::vector<std::vector<size_t> > states_;
};

template<class BP>
class BeliefPropagationProtocolVisitor {
public:
    typedef BP bp_type;
    typedef typename bp_type::value_type value_type;

    BeliefPropagationProtocolVisitor(bool protocolStates = false)
        : protocolStates_(protocolStates)
        { TIC; }
    void operator()(const bp_type& bp)
        {
            double t = TOCN;
            if(runtimes_.size() == 0) { 
                runtimes_.push_back(t);
            }
            else {
                runtimes_.push_back(t + runtimes_[runtimes_.size()-1]);
            }
            
            std::vector<size_t> state;
            bp.arg(state);
            if(protocolStates_) {
                states_.push_back(state);
            }

            value_type value = bp.graphicalModel().evaluate(state);
            values_.push_back(value);

            value_type distance = bp.convergence();
            distances_.push_back(distance);

            std::cout << "step " << runtimes_.size() 
                << ": E=" <<  value
                << ", c=" << distance
                << " after " << runtimes_[runtimes_.size()-1]/1e3 << " s." << std::endl;
            TIC;
        }
    const std::vector<double>& runtimes() const
        { return runtimes_; }
    const std::vector<value_type>& values() const
        { return values_; }
    const std::vector<value_type>& distances() const
        { return distances_; }
    const std::vector<std::vector<size_t> >& states() const
        { return states_; }

private:
    std::vector<double> runtimes_;
    std::vector<value_type> values_;
    std::vector<value_type> distances_;
    bool protocolStates_;
    std::vector<std::vector<size_t> > states_;
};

template<class TRBP>
class TreeReweightedBeliefPropagationProtocolVisitor {
public:
    typedef TRBP trbp_type;
    typedef typename trbp_type::value_type value_type;

    TreeReweightedBeliefPropagationProtocolVisitor(bool protocolStates = false)
        : protocolStates_(protocolStates)
        { TIC; }
    void operator()(const trbp_type& trbp)
        {
            double t = TOCN;
            if(runtimes_.size() == 0) { 
                runtimes_.push_back(t);
            }
            else {
                runtimes_.push_back(t + runtimes_[runtimes_.size()-1]);
            }

            std::vector<size_t> state;
            trbp.arg(state);
            if(protocolStates_) {
                states_.push_back(state);
            }

            value_type value = trbp.graphicalModel().evaluate(state);
            values_.push_back(value);

            value_type distance = trbp.convergence();
            distances_.push_back(distance);

            std::cout << "step " << runtimes_.size() 
                << ": E=" <<  value
                << ", c=" << distance
                << " after " << runtimes_[runtimes_.size()-1]/1e3 << " s." << std::endl;
            TIC;
        }
    const std::vector<double>& runtimes() const
        { return runtimes_; }
    const std::vector<value_type>& values() const
        { return values_; }
    const std::vector<value_type>& distances() const
        { return distances_; }
    const std::vector<std::vector<size_t> >& states() const
        { return states_; }

private:
    std::vector<double> runtimes_;
    std::vector<value_type> values_;
    std::vector<value_type> distances_;
    bool protocolStates_;
    std::vector<std::vector<size_t> > states_;
};

template<class ASTAR, bool Verbose>
class ASTARProtocolVisitor {
public:
    typedef ASTAR astar_type;
    typedef typename astar_type::value_type value_type;

    ASTARProtocolVisitor()
        { TIC; }
    void operator()(
        const astar_type& astar,
        const std::vector<size_t>& conf,
        const size_t& heapsize,
        const value_type& lowerBound, 
        const value_type& upperBound,
        const double& runtime) 
        {
            double t = TOCN;
            if(runtimes_.size() == 0) { 
                runtimes_.push_back(t);
            }
            else {
                runtimes_.push_back(t + runtimes_[runtimes_.size()-1]);
            }
            lowerBounds_.push_back(lowerBound);
            heapSizes_.push_back(heapsize);
            if(conf.size()==astar.graphicalModel().space().dimension())
                upperBounds_.push_back(astar.graphicalModel().evaluate(conf));
            else
                upperBounds_.push_back(upperBound);

            if(Verbose){
                std::cout << "time = " << runtime/1000.0 << "sec"  
                    << "  heapsize = " << heapsize 
                    << "    " <<  lowerBounds_.back() << " <= E(x) <= " << upperBounds_.back()
                    << std::endl;
            }
            TIC; 
        }
    const std::vector<double>& runtimes() const
        { return runtimes_; }
    const std::vector<value_type>& upperBounds() const
        { return upperBounds_; }
    const std::vector<value_type>& lowerBounds() const
        { return lowerBounds_; } const 
    std::vector<size_t>& heapSizes() const
        { return heapSizes_; }

private:
    std::vector<double> runtimes_;
    std::vector<value_type> lowerBounds_;
    std::vector<value_type> upperBounds_;  
    std::vector<size_t> heapSizes_;  
};

struct CommandLineParser {
    enum Algorithm {
        ASTAR, 
        BELIEF_PROPAGATION, 
        GRAPHCUT_PUSH_RELABEL, 
        GRAPHCUT_EDMONDS_KARP,
        GRAPHCUT_KOLMOGOROV,
        ICM, 
        LAZY_FLIPPER, 
        TRBP
    };

    CommandLineParser(int argc, char** argv)
    {
        TCLAP::ValueArg<std::string> outputDatasetArg("D", "output-dataset", "output dataset", false, "state", "string");
        TCLAP::ValueArg<std::string> outputFilenameArg("o", "output-filename", "output file", false, "out.h5", "string");
        TCLAP::ValueArg<std::string> gmDatasetArg("d", "graphical-model-dataset", "graphical model dataset", false, "graphical-model", "string");
        TCLAP::ValueArg<std::string> gmFilenameArg("g", "graphical-model-file", "graphical model file", true, "", "string");

        TCLAP::SwitchArg verboseArg("v", "verbose", "verbose mode");
        TCLAP::SwitchArg protocolArg("p", "protocol", "protocol mode");
        TCLAP::SwitchArg protocolStatesArg("", "protocol-states", "protocol states in protocol mode");

        std::vector<std::string> algorithmAllowedStrings;
        algorithmAllowedStrings.push_back("astar");
        algorithmAllowedStrings.push_back("bp");
        #ifdef WITH_BOOST
            algorithmAllowedStrings.push_back("gc-pr");
            algorithmAllowedStrings.push_back("gc-ek");
            algorithmAllowedStrings.push_back("gc-k");
        #endif
        algorithmAllowedStrings.push_back("icm");
        algorithmAllowedStrings.push_back("lf");
        algorithmAllowedStrings.push_back("trbp");
        TCLAP::ValuesConstraint<std::string> algorithmConstaint(algorithmAllowedStrings);
        TCLAP::ValueArg<std::string> algorithmArg("a", "algorithm", "Algorithm", false, "bp", &algorithmConstaint);

        TCLAP::ValueArg<size_t> bpStepsArg("", "bp-steps", "maximum number of message passing steps", false, 100, "integer");
        TCLAP::ValueArg<float> bpDampingArg("", "bp-damping", "message damping", false, 0.0f, "float");

        TCLAP::ValueArg<size_t> lfDepthArg("", "lf-depth", "maximum size of connected subgraphs.", false, 4, "integer");
        TCLAP::ValueArg<std::string> lfInitFileArg("", "lf-init-file", "file containing the initial state for the Lazy Flipper", false, "", "string");
        TCLAP::ValueArg<std::string> lfInitDatasetArg("", "lf-init-dataset", "file containing the initial state for the Lazy Flipper", false, "state", "string");

        TCLAP::ValueArg<size_t> trbpStepsArg("", "trbp-steps", "maximum number of message passing steps", false, 100, "integer");
        TCLAP::ValueArg<float> trbpDampingArg("", "trbp-damping", "message damping", false, 0.0f, "float");
        TCLAP::ValueArg<std::string> trbpRhoArg("", "trbp-rho", "factor appearence vector", false, "", "string");

        TCLAP::CmdLine cmd("opengm", ' ', "0.1");

        cmd.add(trbpStepsArg);
        cmd.add(trbpDampingArg);
        cmd.add(trbpRhoArg);
        cmd.add(bpStepsArg);
        cmd.add(bpDampingArg);
        cmd.add(lfInitFileArg);
        cmd.add(lfInitDatasetArg);
        cmd.add(lfDepthArg);
        cmd.add(protocolStatesArg);
        cmd.add(protocolArg);
        cmd.add(verboseArg);
        cmd.add(outputDatasetArg);
        cmd.add(outputFilenameArg);
        cmd.add(gmDatasetArg);
        cmd.add(gmFilenameArg);
        cmd.add(algorithmArg);

        cmd.parse(argc, argv);

        outputFileName = outputFilenameArg.getValue();
        outputDatasetName = outputDatasetArg.getValue();
        gmFileName = gmFilenameArg.getValue();
        gmDatasetName = gmDatasetArg.getValue();
        verboseMode = verboseArg.getValue();
        protocolMode = protocolArg.getValue() | protocolStatesArg.getValue();
        protocolStates = protocolStatesArg.getValue();
        if(algorithmArg.getValue().compare("astar") == 0) {
            algorithm = ASTAR;
        }
	else if(algorithmArg.getValue().compare("bp") == 0) {
            algorithm = BELIEF_PROPAGATION;
            bpSteps = bpStepsArg.getValue();
            bpDamping = bpDampingArg.getValue();
        }
        else if(algorithmArg.getValue().compare("gc-pr") == 0) {
            algorithm = GRAPHCUT_PUSH_RELABEL;
        }
        else if(algorithmArg.getValue().compare("gc-ek") == 0) {
            algorithm = GRAPHCUT_EDMONDS_KARP;
        }
        else if(algorithmArg.getValue().compare("gc-k") == 0) {
            algorithm = GRAPHCUT_KOLMOGOROV;
        }
        else if(algorithmArg.getValue().compare("icm") == 0) {
            algorithm = ICM;
        }
        else if(algorithmArg.getValue().compare("lf") == 0) {
            algorithm = LAZY_FLIPPER;
            lfDepth = lfDepthArg.getValue();
            if(lfInitFileArg.isSet()) {
                lfInit = true;
                lfInitFileName = lfInitFileArg.getValue();
                lfInitDatasetName = lfInitDatasetArg.getValue();
            }
            else {
                lfInit = false;
            }
        }
        else if(algorithmArg.getValue().compare("trbp") == 0) {
            algorithm = TRBP;
            trbpSteps = trbpStepsArg.getValue();
            trbpDamping = trbpDampingArg.getValue(); 
	    if(trbpRhoArg.isSet()) {
                rhoFileName = trbpRhoArg.getValue();
		rhoInit     = true;
            }
	    else{
	      rhoInit = false;
	    }
        }
    }

    std::string outputFileName;
    std::string outputDatasetName;
    std::string gmFileName;
    std::string gmDatasetName;
    bool verboseMode;
    bool protocolMode;
    bool protocolStates;
    Algorithm algorithm;
    size_t bpSteps;
    float bpDamping;
    size_t lfDepth;
    bool lfInit;
    std::string lfInitFileName;
    std::string lfInitDatasetName;
    size_t trbpSteps;
    float trbpDamping;
    bool rhoInit;
    std::string rhoFileName;
};

template<class GraphicalModel>
void optimize
(
    const CommandLineParser& clp,
    GraphicalModel& gm
)
{
    typedef typename GraphicalModel::factor_type Factor;
    typedef typename GraphicalModel::value_type value_type;

    // create output file
    hid_t outputFileHandle = marray::hdf5::createFile(clp.outputFileName);

    // optimize
    std::vector<size_t> optimalState(gm.space().dimension());
    if(clp.algorithm == CommandLineParser::BELIEF_PROPAGATION) {
        typedef opengm::BeliefPropagation<GraphicalModel, opengm::Minimizer, opengm::MaxDistance> BP;
        std::cout << "setting up belief propagation... " << std::endl;
        typename BP::Parameter para;
        para.maximumNumberOfSteps_ = clp.bpSteps;
        para.damping_ = clp.bpDamping;
        BP bp(gm, para);
        std::cout << "optimizing... " << std::endl;
        if(clp.protocolMode) {
            BeliefPropagationProtocolVisitor<BP> visitor(clp.protocolStates);
            bp.infer(visitor);
            marray::hdf5::save(outputFileHandle, "runtimes", visitor.runtimes());
            marray::hdf5::save(outputFileHandle, "values", visitor.values());
            marray::hdf5::save(outputFileHandle, "distances", visitor.distances());
            if(clp.protocolStates) {
                marray::Matrix<size_t> states(visitor.states().size(), gm.space().dimension());
                for(size_t j=0; j<visitor.states().size(); ++j) {
                    for(size_t k=0; k<gm.space().dimension(); ++k) {
                        states(j, k) = visitor.states()[j][k];
                    }
                }
                marray::hdf5::save(outputFileHandle, "states", states);
            }
        }
        else if(clp.verboseMode) {
            opengm::BeliefPropagationVerboseVisitor<BP> visitor;
            bp.infer(visitor);
        }
        else {
            bp.infer();
        }
        bp.arg(optimalState);
    }
    else if(clp.algorithm == CommandLineParser::TRBP) {
        typedef opengm::TreeReweightedBeliefPropagation<GraphicalModel, opengm::Minimizer, opengm::MaxDistance> TRBP;
        std::cout << "setting up tree-reweighted belief propagation... " << std::endl;
        typename TRBP::Parameter para;
        para.maximumNumberOfSteps_ = clp.trbpSteps;
        para.damping_ = clp.trbpDamping;
        if(clp.rhoInit) {
            marray::Vector<value_type> rho(gm.numberOfFactors());
            std::cout << "loading edge appearence parameter file... " << std::endl;
            hid_t rhoFile = marray::hdf5::openFile(clp.rhoFileName);
            marray::hdf5::load(rhoFile, "rho", rho);
            marray::hdf5::closeFile(rhoFile);
            // conformity check
            if(rho.size() != gm.numberOfFactors()) {
                throw std::runtime_error("number of entries in rho does not match the number of factors in the graphical model.");
            }
            para.rho_.resize(rho.size());
            for(size_t i=0; i<rho.size(); ++i) {
                para.rho_[i] = rho[i];
            }
        }
        TRBP trbp(gm, para);
        std::cout << "optimizing... " << std::endl;
        if(clp.protocolMode) {
            TreeReweightedBeliefPropagationProtocolVisitor<TRBP> visitor(clp.protocolStates);
            trbp.infer(visitor);
            marray::hdf5::save(outputFileHandle, "runtimes", visitor.runtimes());
            marray::hdf5::save(outputFileHandle, "values", visitor.values());
            marray::hdf5::save(outputFileHandle, "distances", visitor.distances());
            if(clp.protocolStates) {
                marray::Matrix<size_t> states(visitor.states().size(), gm.space().dimension());
                for(size_t j=0; j<visitor.states().size(); ++j) {
                    for(size_t k=0; k<gm.space().dimension(); ++k) {
                        states(j, k) = visitor.states()[j][k];
                    }
                }
                marray::hdf5::save(outputFileHandle, "states", states);
            }
        }
        else if(clp.verboseMode) {
            opengm::TreeReweightedBeliefPropagationVerboseVisitor<TRBP> visitor;
            trbp.infer(visitor);
        }
        else {
            trbp.infer();
        }
        trbp.arg(optimalState);
    }
    else if(clp.algorithm == CommandLineParser::ICM) {
        typedef opengm::ICM<GraphicalModel, opengm::Minimizer> ICM;
        std::cout << "setting up ICM... " << std::endl;
        ICM icm(gm);
        std::cout << "optimizing... " << std::endl;
        if(clp.protocolMode) {
            ICMProtocolVisitor<ICM> visitor(clp.protocolStates);
            icm.infer(visitor);
            marray::hdf5::save(outputFileHandle, "runtimes", visitor.runtimes());
            marray::hdf5::save(outputFileHandle, "values", visitor.values());
            if(clp.protocolStates) {
                marray::Matrix<size_t> states(visitor.states().size(), gm.space().dimension());
                for(size_t j=0; j<visitor.states().size(); ++j) {
                    for(size_t k=0; k<gm.space().dimension(); ++k) {
                        states(j, k) = visitor.states()[j][k];
                    }
                }
                marray::hdf5::save(outputFileHandle, "states", states);
            }
        }
        else if(clp.verboseMode) {
            opengm::ICMVerboseVisitor<ICM> visitor;
            icm.infer(visitor);
        }
        else {
            icm.infer();
        }
        icm.arg(optimalState);
    }
    else if(clp.algorithm == CommandLineParser::LAZY_FLIPPER) {
        typedef opengm::LazyFlipper<GraphicalModel, opengm::Minimizer> LazyFlipper;
        
        // load initial configuration if specified
        marray::Vector<size_t> initialState(gm.space().dimension());
        if(clp.lfInit) {
            std::cout << "loading initial state for the Lazy Flipper... " << std::endl;
            
            hid_t initFile = marray::hdf5::openFile(clp.lfInitFileName);
            marray::hdf5::load(initFile, clp.lfInitDatasetName, initialState);
            marray::hdf5::closeFile(initFile);

            // conformity check
            if(initialState.size() != gm.space().dimension()) {
                throw std::runtime_error("number of entries in initial state does not match the number of variables in the graphical model.");
            }
            for(size_t j=0; j<initialState.size(); ++j) {
                if(initialState[j] >= gm.space().numberOfStates(j)) {
                    throw std::runtime_error("number of states in initial state does not match the number of states of variables in the graphical model.");
                }
            }
        }
        
        std::cout << "setting up the Lazy Flipper... " << std::endl;
        LazyFlipper flipper(gm, clp.lfDepth, initialState.begin());

        std::cout << "optimizing... " << std::endl;
        if(clp.protocolMode) {
            LazyFlipperProtocolVisitor<LazyFlipper> visitor(clp.protocolStates);
            flipper.infer(visitor);
            marray::hdf5::save(outputFileHandle, "runtimes", visitor.runtimes());
            marray::hdf5::save(outputFileHandle, "values", visitor.values());
            marray::hdf5::save(outputFileHandle, "subgraph-sizes", visitor.subgraphSizes());
            marray::hdf5::save(outputFileHandle, "subgraph-forest-sizes", visitor.subgraphForestSizes());
            if(clp.protocolStates) {
                marray::Matrix<size_t> states(visitor.states().size(), gm.space().dimension());
                for(size_t j=0; j<visitor.states().size(); ++j) {
                    for(size_t k=0; k<gm.space().dimension(); ++k) {
                        states(j, k) = visitor.states()[j][k];
                    }
                }
                marray::hdf5::save(outputFileHandle, "states", states);
            }
        }
        else if(clp.verboseMode) {
            opengm::LazyFlipperVerboseVisitor<LazyFlipper> visitor;
            flipper.infer(visitor);
        }
        else {
            flipper.infer();
        }
        flipper.arg(optimalState);
    }
    else if(clp.algorithm == CommandLineParser::ASTAR) {
        typedef opengm::AStar<GraphicalModel, opengm::Minimizer> ASTAR;
        std::cout << "setting up A-Star... " << std::endl;
        typename ASTAR::Parameter para;
        ASTAR astar(gm, para);
        std::cout << "optimizing... " << std::endl;
        if(clp.protocolMode) { 
            if(clp.protocolStates) {
                std::cout << "WARNING: --protocol-states mode not implemented for A-Star search" << std::endl;
            }
            ASTARProtocolVisitor<ASTAR,true> visitor;
            astar.infer(visitor);
            marray::hdf5::save(outputFileHandle, "runtimes",  visitor.runtimes());
            marray::hdf5::save(outputFileHandle, "lower-bounds",    visitor.lowerBounds());
            marray::hdf5::save(outputFileHandle, "upper-bounds", visitor.upperBounds());
            marray::hdf5::save(outputFileHandle, "heapsizes", visitor.heapSizes());
        }
        else if(clp.verboseMode) { 
            opengm::AStarVisitor<ASTAR,true> visitor;
            astar.infer(visitor);
        }
        else {
            astar.infer();
        }

        astar.arg(optimalState);
    }
    #ifdef WITH_BOOST
    else if(clp.algorithm == CommandLineParser::GRAPHCUT_PUSH_RELABEL) {
        typedef opengm::GraphCut<Factor> GraphCut;
        std::cout << "setting up push relabel graph cut... " << std::endl;
        GraphCut gc(gm, GraphCut::PUSH_RELABEL);
        std::cout << "optimizing... " << std::endl;
        std::vector<double> runtimes(1);
        std::vector<value_type> values(1);
        TIC;
        gc.infer();
        runtimes[0] = TOCN;
        gc.arg(optimalState);
        values[0] = gm.evaluate(optimalState.begin());
        marray::hdf5::save(outputFileHandle, "runtimes", runtimes);
        marray::hdf5::save(outputFileHandle, "values", values);
        if(clp.verboseMode || clp.protocolMode) {
            std::cout << "E=" << values[0] << std::endl;
        }
    }
    else if(clp.algorithm == CommandLineParser::GRAPHCUT_EDMONDS_KARP) {
        typedef opengm::GraphCut<Factor> GraphCut;
        std::cout << "setting up Edmonds-Karp graph cut... " << std::endl;
        GraphCut gc(gm, GraphCut::EDMONDS_KARP);
        std::cout << "optimizing... " << std::endl;
        std::vector<double> runtimes(1);
        std::vector<value_type> values(1);
        TIC;
        gc.infer();
        runtimes[0] = TOCN;
        gc.arg(optimalState);
        values[0] = gm.evaluate(optimalState.begin());
        marray::hdf5::save(outputFileHandle, "runtimes", runtimes);
        marray::hdf5::save(outputFileHandle, "values", values);
        if(clp.verboseMode || clp.protocolMode) {
            std::cout << "E=" << values[0] << std::endl;
        }
    }
    else if(clp.algorithm == CommandLineParser::GRAPHCUT_KOLMOGOROV) {
        typedef opengm::GraphCut<Factor> GraphCut;
        std::cout << "setting up Kolmogorov graph cut... " << std::endl;
        GraphCut gc(gm, GraphCut::KOLMOGOROV);
        std::cout << "optimizing... " << std::endl;
        std::vector<double> runtimes(1);
        std::vector<value_type> values(1);
        TIC;
        gc.infer();
        runtimes[0] = TOCN;
        gc.arg(optimalState);
        values[0] = gm.evaluate(optimalState.begin());
        marray::hdf5::save(outputFileHandle, "runtimes", runtimes);
        marray::hdf5::save(outputFileHandle, "values", values);
        if(clp.verboseMode || clp.protocolMode) {
            std::cout << "E=" << values[0] << std::endl;
        }
    }
    #endif

    std::cout << "saving output... " << std::endl;
    marray::hdf5::save(outputFileHandle, clp.outputDatasetName, optimalState);
    H5Fclose(outputFileHandle); 
    /*
	std::cout <<"E[ ";
	for(size_t i=0; i<optimalState.size();++i)
	  std::cout << optimalState[i] << " ";
	std::cout << "] = "<< gm.evaluate( optimalState) << std::endl;
    */

}

// main

int main(int argc, char** argv) {
    std::cout << std::endl
        << "OpenGM Optimizer v1.0. Copyright (c) 2010 by B. Andres & J. H. Kappes" << std::endl 
        << "bjoern.andres@iwr.uni-heidelberg.de, kappes@math.uni-heidelberg.de" << std::endl 
        << "HCI, IWR, University of Heidelberg." << std::endl 
        << std::endl
        << "Please cite: " << std::endl
        << "B. Andres, J. H. Kappes, U. Koethe, C. Schnoerr, F. A. Hamprecht." << std::endl
        << "An Empirical Comparison of Inference Algorithms for Graphical Models" << std::endl
        << "with Higher Order Factors Using OpenGM. DAGM 2010." << std::endl
        << std::endl;

    try {
        CommandLineParser clp(argc, argv);

        // load serialized graphical model
        typedef double value_type; 
        std::cout << "loading model... " << std::endl;
        marray::Vector<value_type> serializedModel;
        {
            hid_t gmFileHandle = marray::hdf5::openFile(clp.gmFileName);
            marray::hdf5::load(gmFileHandle, clp.gmDatasetName, serializedModel);
            marray::hdf5::closeFile(gmFileHandle); 
        }
        
        // deserialize and optimize model
        bool binaryVariablesOnly = true;
        size_t numberOfVariables = static_cast<size_t>(serializedModel[0]);
        for(size_t j=2; j<numberOfVariables+1; ++j) {
            if(serializedModel[j] != 2) {
                binaryVariablesOnly = false;
                break;
            }
        }
        if(binaryVariablesOnly) {
            typedef opengm::DiscreteSpace Space;
            typedef opengm::ExplicitFactor<value_type, ArrayPow2<value_type> > Factor;
            typedef opengm::GraphicalModel<Factor, opengm::Adder> GraphicalModel;

            if(clp.verboseMode || clp.protocolMode) {
                std::cout << "switching to high-performance mode for binary variables." << std::endl;
            }
            Space space;
            GraphicalModel gm;
            opengm::deserialize(serializedModel, gm, space);

            optimize<GraphicalModel>(clp, gm);
        }
        else {
            typedef opengm::DiscreteSpace Space;
            typedef opengm::ExplicitFactor<value_type> Factor;
            typedef opengm::GraphicalModel<Factor, opengm::Adder> GraphicalModel;

            Space space;
            GraphicalModel gm;
            opengm::deserialize(serializedModel, gm, space);

            optimize<GraphicalModel>(clp, gm);
        }
        return 0;
    }
    catch(TCLAP::ArgException& e) {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl; 
    }
    catch(std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
