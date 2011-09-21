/// UFD. Copyright (c) 2010 by Bjoern Andres
///
/// This software was developed by Bjoern Andres 
/// Enquiries shall be directed to: bjoern.andres@iwr.uni-heidelberg.de
///
/// All advertising materials mentioning features or use of this software must
/// display the following acknowledgement: ``This product includes UFD developed 
/// by Bjoern Andres. Please direct enquiries concerning UFD to 
/// bjoern.andres@iwr.uni-heidelberg.de''.
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
///   display the following acknowledgement: ``This product includes UFD developed 
///   by Bjoern Andres. Please direct enquiries concerning UFD to 
///   bjoern.andres@iwr.uni-heidelberg.de''.
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
#ifndef UFD_HXX
#define UFD_HXX

#include <cassert>
#include <vector>
#include <map>

namespace ufd {

template<class T = size_t>
class Partition {
public:
	Partition();
	Partition(const T&);

	// query
	T find(T);
	inline T numberOfElements() const {return numberOfElements_;}
	inline T numberOfSets() const {return numberOfSets_;}
	template<class Iterator>
		inline void elementLabeling(Iterator);
	template<class Iterator>
		inline void representatives(Iterator);
	inline void representativeLabeling(std::map<T, T>&);

	// manipulation
	void merge(T, T);
	void insert(const T&);

private:
	std::vector<T> parents_;
	std::vector<T> ranks_;
	T numberOfElements_;
    T numberOfSets_;
};

template<class T>
Partition<T>::Partition()
: parents_(std::vector<T>()),
  ranks_(std::vector<T>()),
  numberOfElements_(0),
  numberOfSets_(0)
{}

/**
 * \param size 
 */ 
template<class T>
Partition<T>::Partition
(
	const T& size
)
: parents_(std::vector<T>(size_t(size))),
  ranks_(std::vector<T>(size_t(size))),
  numberOfElements_(size),
  numberOfSets_(size)
{
	for(T j=0; j<size; ++j) {
		parents_[size_t(j)] = j;
	}
}

template<class T>
T
Partition<T>::find
(
	T element // copy to work with
)
{
	assert(element < numberOfElements());
	// find the root
	T root = element;
	while(parents_[size_t(root)] != root) {
		root = parents_[size_t(root)];
	}
	// path compression
	while(element != root) {
		T tmp = parents_[size_t(element)];
		parents_[size_t(element)] = root;
		element = tmp;
	}
	return root;
}

template<class T>
void Partition<T>::merge
(
	T element1,
	T element2
)
{
	assert(element1 < numberOfElements());
	assert(element2 < numberOfElements());
	// merge by rank
	element1 = find(element1);
	element2 = find(element2);
	if(ranks_[size_t(element1)] < ranks_[size_t(element2)]) {
		parents_[size_t(element1)] = element2;
		--numberOfSets_;
	}
	else if(ranks_[size_t(element1)] > ranks_[size_t(element2)]) {
		parents_[size_t(element2)] = element1;
		--numberOfSets_;
	}
	else if(element1 != element2) {
		parents_[size_t(element2)] = element1;
		++ranks_[size_t(element1)];
		--numberOfSets_;
	}
}

/** insert new sets, each consisting of one new element
*/
template<class T>
void Partition<T>::insert
(
	const T& number
)
{
	assert(number >= 0);

	ranks_.insert(ranks_.end(), size_t(number), T(0));
	parents_.insert(parents_.end(), size_t(number), T(0));
	for(T j=numberOfElements_; j<numberOfElements_+number; ++j) {
		parents_[size_t(j)] = j;
	}
	numberOfElements_ += number;
	numberOfSets_ += number;
}

/** Representatives

Output all elements which are set representatives
to a container through the given iterator
*/
template<class T>
template<class Iterator>
void Partition<T>::representatives
(
	Iterator it
)
{
	for(T j=0; j<numberOfElements(); ++j) {
		if(parents_[size_t(j)] == j) {
			*it = j;
			++it;
		}
	}
}

/** Dense labeling of set representatives
 *  
 */
template<class T>
void Partition<T>::representativeLabeling
(
	std::map<T, T>& out
)
{
	out.clear();	
	std::vector<T> r( (size_t) (numberOfSets()) ); 
	representatives(r.begin());
	for(T j=0; j<numberOfSets(); ++j) {
		out[ r[size_t(j)] ] = j;
	}
}

/** dense labeling of all elements
*/
template<class T>
template<class Iterator>
void Partition<T>::elementLabeling
(
	Iterator out
)
{
	std::map<T, T> rl;
	representativeLabeling(rl);
	for(T j=0; j<numberOfElements(); ++j) {
		*out = rl[find(j)];
		++out;
	}
}

} // namespace ufd

#endif //UFD_HXX
