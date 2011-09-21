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
#ifndef OPENGM_GRAPHCUT_HXX
#define OPENGM_GRAPHCUT_HXX

#include <queue>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/kolmogorov_max_flow.hpp>

#include "opengm/explicitfactor.hxx"
#include "opengm/graphicalmodel.hxx"
#include "opengm/adder.hxx"
#include "opengm/minimizer.hxx"
#include "opengm/inference/inference.hxx"

namespace opengm {

// Graph Cut Optimizer.
//
// This optimizer minimizes the value of a graphical model
// of the type GraphicalModel<Factor, opengm::Adder>
// 
template<class Factor>
class GraphCut 
: Inference<GraphicalModel<Factor, opengm::Adder>, 
            opengm::Minimizer> 
{
public:
    typedef Factor factor_type;
    typedef typename factor_type::space_type space_type; 
    typedef typename factor_type::value_type value_type;
    typedef typename space_type::state_type state_type;
    typedef GraphicalModel<Factor, opengm::Adder> gm_type;
    enum MaxFlowAlgorithm { PUSH_RELABEL, EDMONDS_KARP, KOLMOGOROV };

    // construction
    GraphCut(const gm_type&, const MaxFlowAlgorithm& = PUSH_RELABEL);

    // query
    std::string name() const;
    const gm_type& graphicalModel() const;

    // manipulation
    void setMaxFlowAlgorithm(const MaxFlowAlgorithm&);

    // inference
    InferenceTermination infer();
    InferenceTermination arg(std::vector<state_type>&, const size_t& = 1) const;

private:
    // boost graph library
    typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> graph_traits;
    typedef graph_traits::edge_descriptor edge_descriptor;
    typedef graph_traits::vertex_descriptor vertex_descriptor;
    struct Edge {
        Edge() : capacity(value_type()), residual(value_type()), reverse(edge_descriptor()) {}
        value_type capacity;
        value_type residual;
        edge_descriptor reverse;
    };
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, size_t, Edge> graph_type;
    typedef typename boost::graph_traits<graph_type>::edge_iterator edge_iterator;
    typedef typename boost::graph_traits<graph_type>::out_edge_iterator out_edge_iterator;

    void constructGraph();
    void printGraph();
    void cutPushRelabel();
    void cutEdmondsKarp();
    void cutKolmogorov();

    const gm_type& gm_;
    MaxFlowAlgorithm maxFlowAlgorithm_;
    graph_type graph_;
    std::vector<state_type> state_;
};

// public interface

template<class Factor>
GraphCut<Factor>::GraphCut
(
    const typename GraphCut::gm_type& gm,
    const MaxFlowAlgorithm& maxFlowAlgorithm
)
: gm_(gm),
  maxFlowAlgorithm_(maxFlowAlgorithm)
{
    for(size_t j=0; j<gm_.space().dimension(); ++j) {
        if(gm_.space().numberOfStates(j) != 2) {
            throw std::runtime_error("This implementation of the graph cut optimizer supports only binary variables.");
        }
    }
    for(size_t j=0; j<gm_.numberOfFactors(); ++j) {
        if(gm_[j].numberOfVariables() > 2) {
            throw std::runtime_error("This implementation of the graph cut optimizer supports only factors of order <= 2.");
        }
        if(!gm_[j].isSubmodular()) {
            throw std::runtime_error("This implementation of the graph cut optimizer supports only submodular factors.");
        }
    }
}

template<class Factor>
inline std::string 
GraphCut<Factor>::name() const
{
    return "GraphCut";
}

template<class Factor>
inline const typename GraphCut<Factor>::gm_type& 
GraphCut<Factor>::graphicalModel() const
{
    return gm_;
}

template<class Factor>
inline void 
GraphCut<Factor>::setMaxFlowAlgorithm
(
    const MaxFlowAlgorithm& maxFlowAlgorithm
)
{
    maxFlowAlgorithm_ = maxFlowAlgorithm;
}

template<class Factor>
inline InferenceTermination 
GraphCut<Factor>::infer()
{
    constructGraph();
    if(maxFlowAlgorithm_ == PUSH_RELABEL) {
        cutPushRelabel();
    }
    else if(maxFlowAlgorithm_ == EDMONDS_KARP) {
        cutEdmondsKarp();
    }
    else if(maxFlowAlgorithm_ == KOLMOGOROV) {
        cutKolmogorov();
    }
    return NORMAL;
}

template<class Factor>
inline InferenceTermination 
GraphCut<Factor>::arg
(
    std::vector<state_type>& arg, 
    const size_t& n
) const
{
    if(n > 1) {
        return UNKNOWN;
    }
    else {
        // skip source and sink
        arg.resize(state_.size() - 2);
        for(size_t j=2; j<state_.size(); ++j) {
            arg[j-2] = state_[j];
        }
        return NORMAL;
    }
}

// private member functions

template<class Factor>
void
GraphCut<Factor>::constructGraph()
{
    graph_ = graph_type(gm_.space().dimension()+2); 
    // add source and sink vertex
    for(size_t j=0; j<gm_.numberOfFactors(); ++j) {
        const factor_type& factor = gm_[j];
        if(factor.numberOfVariables() == 0) {
            // constant factors
            // ignoring
        }
        else if(factor.numberOfVariables() == 1) {
            // 1st order factor
            if(factor(0) <= factor(1)) {
                // edge
                std::pair<edge_descriptor, bool> e =
                    add_edge(0, factor.variableIndex(0)+2, graph_); // from source
                graph_[e.first].capacity = factor(1) - factor(0); // cost
                // reverse edge
                std::pair<edge_descriptor, bool> er =
                    add_edge(factor.variableIndex(0)+2, 0, graph_); 
                graph_[e.first].reverse = er.first;
                graph_[er.first].reverse = e.first; 
            }
            else {
                // edge
                std::pair<edge_descriptor, bool> e =
                    add_edge(factor.variableIndex(0)+2, 1, graph_); // to sink
                graph_[e.first].capacity = factor(0) - factor(1); // cost
                // reverse edge
                std::pair<edge_descriptor, bool> er =
                    add_edge(1, factor.variableIndex(0)+2, graph_); 
                graph_[e.first].reverse = er.first;
                graph_[er.first].reverse = e.first; 
            }
        }
        else if(factor.numberOfVariables() == 2) {
            // 2nd order factor
            const value_type& A = factor(0,0);
            const value_type& B = factor(0,1);
            const value_type& C = factor(1,0);
            const value_type& D = factor(1,1);
            // first variabe
            if(C > A) {
                // edge
                std::pair<edge_descriptor, bool> e =
                    add_edge(0, factor.variableIndex(0)+2, graph_); // from source
                graph_[e.first].capacity = C - A; // cost
                // reverse edge
                std::pair<edge_descriptor, bool> er =
                    add_edge(factor.variableIndex(0)+2, 0, graph_); 
                graph_[e.first].reverse = er.first;
                graph_[er.first].reverse = e.first; 
            }
            else if(C < A) {
                // edge
                std::pair<edge_descriptor, bool> e = 
                    add_edge(factor.variableIndex(0)+2, 1, graph_); // to sink
                graph_[e.first].capacity = A - C; // cost
                // reverse edge
                std::pair<edge_descriptor, bool> er = 
                    add_edge(1, factor.variableIndex(0)+2, graph_); 
                graph_[e.first].reverse = er.first;
                graph_[er.first].reverse = e.first; 
            }
            // second variable
            if(D > C) {
                // edge
                std::pair<edge_descriptor, bool> e = 
                    add_edge(0, factor.variableIndex(1)+2, graph_); // from source
                graph_[e.first].capacity = D - C; // cost
                // reverse edge
                std::pair<edge_descriptor, bool> er = 
                    add_edge(factor.variableIndex(1)+2, 0, graph_); 
                graph_[e.first].reverse = er.first;
                graph_[er.first].reverse = e.first; 
            }
            else if(D < C) {
                // edge
                std::pair<edge_descriptor, bool> e = 
                    add_edge(factor.variableIndex(1)+2, 1, graph_); // to sink
                graph_[e.first].capacity = C - D; // cost
                // reverse edge
                std::pair<edge_descriptor, bool> er = 
                    add_edge(1, factor.variableIndex(1)+2, graph_); 
                graph_[e.first].reverse = er.first;
                graph_[er.first].reverse = e.first; 
            }
            // submodular term
            if(B + C - A - D > 0) {
                // edge
                std::pair<edge_descriptor, bool> e = 
                    add_edge(factor.variableIndex(0)+2, factor.variableIndex(1)+2, graph_); 
                graph_[e.first].capacity = B + C - A - D; // cost
                // reverse edge
                std::pair<edge_descriptor, bool> er = 
                    add_edge(factor.variableIndex(1)+2, factor.variableIndex(0)+2, graph_); 
                graph_[e.first].reverse = er.first;
                graph_[er.first].reverse = e.first; 
            }
        }
        else {
            // higher order factor
            throw std::runtime_error("This implementation of the graph cut optimizer does not support factors of order >2.");
        }
    }
}

template<class Factor>
void GraphCut<Factor>::cutEdmondsKarp()
{
    // compute max flow
    std::vector<boost::default_color_type> color(num_vertices(graph_));
    std::vector<edge_descriptor> pred(num_vertices(graph_));
    /*T flow = */edmonds_karp_max_flow(graph_, 0, 1,
        get(&Edge::capacity, graph_),
        get(&Edge::residual, graph_),
        get(&Edge::reverse, graph_),
        &color[0], &pred[0]
    );

    // find (s,t)-cut set
    state_ = std::vector<state_type>(num_vertices(graph_));
    for(size_t j=2; j<num_vertices(graph_); ++j) {
        if(color[j] == boost::black_color) {
            state_[j] = 0;
        }
        else if(color[j] == boost::white_color) {
            state_[j] = 1;
        }
        else {
            throw std::runtime_error("At least one vertex is labeled neither black nor white.");
        }
    }
}

template<class Factor>
void GraphCut<Factor>::cutPushRelabel()
{
    // compute max flow
    /*T flow = */push_relabel_max_flow(graph_, 0, 1, 
        get(&Edge::capacity, graph_), get(&Edge::residual, graph_),
        get(&Edge::reverse, graph_), get(boost::vertex_index_t(), graph_));

    // find (s,t)-cut set
    state_ = std::vector<state_type>(num_vertices(graph_), 1);
    state_[0] = 0; // source
    state_[1] = 0; // sink
    typedef typename boost::property_map<graph_type, boost::vertex_index_t>::type VertexIndexMap;
    VertexIndexMap vertexIndexMap = get(boost::vertex_index, graph_);
    std::queue<vertex_descriptor> q;
    q.push(*(vertices(graph_).first)); // source
    while(!q.empty()) {
        out_edge_iterator current, end;
        tie(current, end) = out_edges(q.front(), graph_); 
        q.pop();
        while(current != end) {
            if(graph_[*current].residual > 0) {
                vertex_descriptor v = target(*current, graph_);
                if(vertexIndexMap[v] > 1 && state_[vertexIndexMap[v]] == 1) {
                    state_[vertexIndexMap[v]] = 0;
                    q.push(v);
                }
            }
            ++current;
        }
    }
}

template<class Factor>
void GraphCut<Factor>::cutKolmogorov() {
    // compute max flow
    std::vector<boost::default_color_type> color(num_vertices(graph_));
    std::vector<edge_descriptor> pred(num_vertices(graph_));
    std::vector<vertex_descriptor> dist(num_vertices(graph_));
    /*T flow = */kolmogorov_max_flow(graph_,
        get(&Edge::capacity, graph_),
        get(&Edge::residual, graph_),
        get(&Edge::reverse, graph_),
        &pred[0],
        &color[0],
        &dist[0],
        get(boost::vertex_index, graph_),
        0, 1
        );

    // find (s,t)-cut set
    state_ = std::vector<state_type>(num_vertices(graph_));
    for(size_t j=2; j<num_vertices(graph_); ++j) {
        if(color[j] == boost::black_color || color[j] == boost::gray_color) {
            state_[j] = 0;
        }
        else if(color[j] == boost::white_color) {
            state_[j] = 1;
        }
    }
}

template<class Factor>
void
GraphCut<Factor>::printGraph()
{
    edge_iterator it, end;
    for(tie(it, end) = edges(graph_); it != end; ++it) {
        std::cout << '(' << source(*it, graph_) << ", " << target(*it, graph_) << "): "
            << "capacity=" << graph_[*it].capacity << ", residual=" << graph_[*it].residual << '.' << std::endl;
    }
}

} // namespace opengm

#endif // #ifndef OPENGM_GRAPHCUT_HXX
