/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
 *
 *  Copyright 2013 Benjamin Tams
 *
 *  THIMBLE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  THIMBLE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with THIMBLE. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file BinomialIterator.cpp
 *
 * @brief
 *            Implements functionalities provided by 'BinomialIterator.h'
 *            which provides a mechanism for iterating through all choices of
 *            <i>k</i> elements from a vector of size <i>n</i>.
 *
 * @author Benjamin Tams
 */

#include <cstdlib>
#include <iostream>

#include <thimble/math/BinomialIterator.h>

using namespace std;

namespace thimble {

	/**
	* @brief
	*            Creates an instance of an iterator to iterate through
	*            all choices of <i>k</i> entries (of different positions)
	*            from an array of size <i>n</i>.
	*
	* @details
	*            see 'BinomialIterator.h'
	*/
	BaseBinomialIterator::BaseBinomialIterator( int n , int k ) {

		// Check if arguments are reasonable.
		if ( n < 0 || k < 0 || n < k ) {
			cerr << "BinomialIterator: Bad arguments." << endl;
			exit(EXIT_FAILURE);
		}

		// Adopt 'n' and 'k'.
		this->n = n;
		this->k = k;

		// Initially, the least index position from where to choose
		// the first entry from an array is '0'.
		this->idx = 0;

		// If it is possible to choose 'k-1' distinct indices which are
		// at higher position than '0', ...
		if ( n-1 >= k-1 && k-1 > 0 ) {
			// ..., recursively create an iterator to choose 'k-1' from 'n-1'
			// at initial state; ...
			this->subIteratorPtr = new BaseBinomialIterator(n-1,k-1);
		} else {
			// ...; otherwise, there is no sub-iterator.
			this->subIteratorPtr = NULL;
		}
	}

	/**
	 * @brief
	 *            Destructor.
	 *
	 * @details
	 *            Frees all the memory that is helt by the iterator.
	 */
	BaseBinomialIterator::~BaseBinomialIterator() {
		// Recursively, delete the sub-iterators.
		delete this->subIteratorPtr;
	}

	/**
	 * @brief
	 *            Changes the iterator's state to choose the next <i>k</i>
	 *            entries of different positions from an array of size
	 *            <i>n</i>.
	 *
	 * @details
	 *            see 'BinomialIterator.h'
	 */
	bool BaseBinomialIterator::next() {

		// If the sub-iterator is NULL ...
		if ( this->subIteratorPtr == NULL ) {

			// ... we first check whether the iterator can change its
			// state to a valid next state; ...
			if ( this->idx+1 >= this->n || this->k == 0 ) {
				// ...; if not, return 'false' to indicate that the iterator
				// already is at its final state; ...
				return false;
			}

			// ...; otherwise, choose the next least index position...
			this->idx += 1;

			// ... and return 'true' to indicate the iterator successfully
			// changed its state.
			return true;
		}

		// If the sub-iterator is not NULL, we first attempt to recursively
		// change the sub-iterator's state....
		if ( this->subIteratorPtr->next() ) {
			// ... and if successful, return 'true' to indicate that the
			// iterator's state changed successfully; ...
			return true;
		}

		// ...; otherwise, if the sub-iterator already has reached its final
		// state, increase the least index positions that is choosed by the
		// iterator.
		this->idx += 1;

		// If the index is not too large, i.e. if it possible to choose
		// further 'k-1' distinct indices, ...
		if ( this->n-this->idx-1 >= this->k-1 && this->k-1 >= 0 ) {
			// ..., we create a sub-iterator to choose 'k-1' out of the
			// remaining  n-idx-1' index positions larger than 'idx'.
			delete this->subIteratorPtr;
			this->subIteratorPtr = new BaseBinomialIterator
					(this->n-this->idx-1,this->k-1);
		} else {
			// Otherwise, if the index is too large, the iterator was already
			// at its final state.
			this->idx -= 1;
			return false;
		}

		return this->idx < this->n;
	}

}



