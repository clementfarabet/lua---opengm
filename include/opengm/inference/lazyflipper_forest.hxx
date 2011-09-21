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
#ifndef OPENGM_LAZYFLIPPER_SUBGRAPHFOREST_HXX
#define OPENGM_LAZYFLIPPER_SUBGRAPHFOREST_HXX

#include <string>

namespace opengm {

/// Forest with Level Order Traversal.
///
/// - no manipulation after construction.
/// - level Successors must be set manually
/// - implementation nor const correct
///
template<class T>
class Forest {
public:
    typedef T Value;
    typedef size_t NodeIndex;
    // static attribute
    static const NodeIndex NONODE = -1;
    typedef size_t Level;

    // constructor
    Forest();

    // forest query
    size_t size();
    size_t levels();
    NodeIndex levelAnchor(const Level&);

    // forest manipulation
    NodeIndex push_back(const Value&, NodeIndex);

    // forest testing and output
    size_t testInvariant();
    std::string asString();

    // node query
    Value& value(NodeIndex);
    Level level(NodeIndex);
    NodeIndex parent(NodeIndex);
    NodeIndex levelOrderSuccessor(NodeIndex);
    size_t numberOfChildren(NodeIndex); 
    NodeIndex child(NodeIndex, const size_t&);

    // node manipulation
    void setLevelOrderSuccessor(NodeIndex, NodeIndex);

private:
    struct Node {
        Node(const Value& value)
            : value_(value), parent_(NONODE), 
            children_(std::vector<NodeIndex>()),
            level_(0), levelOrderSuccessor_(NONODE)
        {}
        Value value_;
        NodeIndex parent_;
        std::vector<NodeIndex> children_;
        Level level_;
        NodeIndex levelOrderSuccessor_;
    };

    std::vector<Node> nodes_;
    std::vector<NodeIndex> levelAnchors_;
};

// implementation

template<class T>
inline Forest<T>::Forest()
:   nodes_(std::vector<typename Forest<T>::Node>()),
    levelAnchors_(std::vector<typename Forest<T>::NodeIndex>())
{
}

template<class T>
inline size_t 
Forest<T>::levels()
{
    return levelAnchors_.size();
}

template<class T>
inline size_t 
Forest<T>::size()
{
    return nodes_.size();
}

template<class T>
inline typename Forest<T>::NodeIndex 
Forest<T>::levelAnchor
(
    const typename Forest<T>::Level& level
)
{
    Assert(NO_DEBUG || level < levels());
    return levelAnchors_[level];
}

template<class T>
inline typename Forest<T>::Value& 
Forest<T>::value
(
    typename Forest<T>::NodeIndex n
)
{
    Assert(NO_ARG_TEST || n < nodes_.size());
    return nodes_[n].value_;
}

template<class T>
inline typename Forest<T>::Level 
Forest<T>::level
(
    typename Forest<T>::NodeIndex n
) 
{
    Assert(NO_ARG_TEST || n < nodes_.size());
    return nodes_[n].level_;
}

template<class T>
inline typename Forest<T>::NodeIndex 
Forest<T>::parent
(
    typename Forest<T>::NodeIndex n
) 
{
    Assert(NO_ARG_TEST || n < nodes_.size());
    return nodes_[n].parent_;
}

template<class T>
inline typename Forest<T>::NodeIndex 
Forest<T>::levelOrderSuccessor
(
    typename Forest<T>::NodeIndex n
) 
{
    Assert(NO_ARG_TEST || n < nodes_.size());
    return nodes_[n].levelOrderSuccessor_;
}

template<class T>
inline size_t 
Forest<T>::numberOfChildren
(
    typename Forest<T>::NodeIndex n
) 
{
    Assert(NO_ARG_TEST || n < nodes_.size());
    return nodes_[n].children_.size();
}

template<class T>
inline typename Forest<T>::NodeIndex 
Forest<T>::child
(
    typename Forest<T>::NodeIndex n,
    const size_t& j
) 
{
    Assert(NO_ARG_TEST || (n<nodes_.size() && j<nodes_[n].children_.size()));
    return nodes_[n].children_[j];
}

template<class T>
typename Forest<T>::NodeIndex 
Forest<T>::push_back
(
    const Value& value,
    typename Forest<T>::NodeIndex parentNodeIndex
)
{
    Assert(NO_ARG_TEST || (parentNodeIndex == NONODE || parentNodeIndex < nodes_.size()));
    // lock here in parallel code
    NodeIndex nodeIndex = nodes_.size();
    {
        Node node(value);
        nodes_.push_back(node);
        // unlock here in parallel code
        Assert(nodes_.size() == nodeIndex + 1); // could fail in parallel code
    }
    if(parentNodeIndex != NONODE) {
        nodes_[nodeIndex].parent_ = parentNodeIndex;
        nodes_[parentNodeIndex].children_.push_back(nodeIndex);
        nodes_[nodeIndex].level_ = nodes_[parentNodeIndex].level_ + 1;
    }
    if(nodes_[nodeIndex].level_ >= levelAnchors_.size()) {
        Assert(NO_DEBUG || levelAnchors_.size() == nodes_[nodeIndex].level_);
        levelAnchors_.push_back(nodeIndex);
    }
    return nodeIndex;
}

// returns the number of root nodes
template<class T>
size_t 
Forest<T>::testInvariant() 
{
    if(nodes_.size() == 0) { 
        // tree is empty
        Assert(levelAnchors_.size() == 0);
        return 0;
    }
    else { 
        // tree is not empty
        Assert(levelAnchors_.size() != 0);
        size_t numberOfRoots = 0;
        size_t nodesVisited = 0;
        Level level = 0;
        NodeIndex p = levelAnchors_[0];
        while(p != NONODE) {
            ++nodesVisited;
            Assert(this->level(p) == level);
            if(level == 0) { 
                // p is a root node index
                Assert(parent(p) == NONODE);
                ++numberOfRoots;
            }
            else { 
                // p is not a root node index
                Assert(parent(p) != NONODE);
                // test if p is among the children of its parent:
                bool foundP = false;
                for(size_t j=0; j<nodes_[parent(p)].children_.size(); ++j) {
                    if(nodes_[parent(p)].children_[j] == p) {
                        foundP = true;
                        break;
                    }
                }
                Assert(foundP);
            }
            // continue traversal in level-order
            if(levelOrderSuccessor(p) != NONODE) {
                p = levelOrderSuccessor(p);
            }
            else {   
                if(level+1 < levelAnchors_.size()) { 
                    // tree has more levels
                    ++level;
                    p = levelAnchors_[level];
                }
                else { 
                    // tree has no more levels
                    break; 
                }
            }
        }
        Assert(nodesVisited == nodes_.size()); 
        Assert(levels() == level+1);
        return numberOfRoots;
    }
}

template<class T>
std::string 
Forest<T>::asString() 
{
    std::ostringstream out(std::ostringstream::out);
    for(size_t level=0; level<levels(); ++level) {
        NodeIndex p = levelAnchor(level);
        while(p != NONODE) {
            // print all variable indices on the path to the root
            NodeIndex q = p;
            while(q != NONODE) {
                // out << value(q) << ' ';
                out << value(q)+1 << ' '; // ??? replace by previous line!!!
                q = parent(q);
            }
            out << std::endl;
            // proceed
            p = levelOrderSuccessor(p);
        }
    }
    return out.str();
}

template<class T>
inline void 
Forest<T>::setLevelOrderSuccessor
(
    typename Forest<T>::NodeIndex nodeIndex, 
    typename Forest<T>::NodeIndex successorNodeIndex
)
{
    Assert(NO_ARG_TEST || (nodeIndex < nodes_.size() && 
        successorNodeIndex < nodes_.size()));
    nodes_[nodeIndex].levelOrderSuccessor_ = successorNodeIndex;
}

} // namespace opengm 

#endif // OPENGM_LAZYFLIPPER_SUBGRAPHFOREST_HXX
