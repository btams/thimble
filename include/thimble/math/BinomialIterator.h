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
 * @file BinomialIterator.h
 *
 * @brief
 *            Provides a mechanism for iterating through all choices of
 *            <i>k</i> elements from a vector of size <i>n</i>.
 *
 * @author Benjamin Tams
 *
 * @see thimble::BinomialIterator
 */

#ifndef THIMBLE_BINOMIALITERATOR_H_
#define THIMBLE_BINOMIALITERATOR_H_

#include <cstdlib>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Base iterator which only provides a mechanism for changing
	 *            the iterator's state but not for choosing <i>k</i> from
	 *            <i>n</i> entries.
	 *
	 * @details
	 *            The template class \link BinomialIterator\endlink extends this
	 *            class by providing an additional function for selecting
	 *            <i>k</i> entries of different positions from an array of size
	 *            <i>n</i>.
	 *
	 * @attention
	 *            Do not use this class for iteration. Use
	 *            \link BinomialIterator\endlink.
	 *
	 * @see BinomialIterator
	 */
	class THIMBLE_DLL BaseBinomialIterator {

	private:

		/**
		 * @brief
		 *            The size of the arrays from where <i>k</i> entries of
		 *            different positions are chosen in each iteration.
		 */
		int n;

		/**
		 * @brief
		 *            The number of entries of different positions that are
		 *            chosen from an array of size <i>n</i> in each iteration.
		 */
		int k;

		/**
		 * @brief
		 *            Always holds the least index position from where to
		 *            choose the first entry.
		 */
		int idx;

		/**
		 * @brief
		 *            Sub-iterator used to define the current state of
		 *            the iterator.
		 *
		 * @details
		 *            If this sub-iterator is not <code>NULL</code>, the
		 *            next least index position from where to choose the
		 *            next entry is
		 *            <pre>
		 *            this->idx+subIteratorPtr->idx+1;
		 *            </pre>
		 *            In such a way, the implementation of our <i>n</i>
		 *            choose <i>k</i> iterator is recursive.
		 */
		BaseBinomialIterator *subIteratorPtr;

	protected:

		/**
		 * @brief
		 *            Creates an instance of an iterator to iterate through
		 *            all choices of <i>k</i> entries (of different positions)
		 *            from an array of size <i>n</i>.
		 *
		 * @param n
		 *            Size of the arrays from where <i>k</i> elements are
		 *            chosen in each iteration.
		 *
		 * @param k
		 *            Specifies the number of elements that are chosen from
		 *            an array of size <i>n</i> in each iteration.
		 *
		 * @warning
		 *            If not sufficient memory could be provided, exceptions
		 *            may be thrown correspondingly.
		 *
		 * @warning
		 *            If <i>n</i> or <i>k</i> is smaller than zero or if
		 *            <i>n</i> is smaller than <i>k</i> an error message will
		 *            be written to <code>stderr</code> and an exit with
		 *            status 'EXIT_FAILURE' will be caused.
		 */
		BaseBinomialIterator( int n , int k );

	public:

		/**
		 * @brief
		 *            Destructor.
		 *
		 * @details
		 *            Frees all the memory that is helt by the iterator.
		 */
		~BaseBinomialIterator();

		/**
		 * @brief
		 *            Changes the iterator's state to choose the next <i>k</i>
		 *            entries of different positions from an array of size
		 *            <i>n</i>.
		 *
		 * @details
		 *            If the state of the iterator successfully changed, the
		 *            result of the function will be <code>true</code>.
		 *            Otherwise, if the iterator has reached its final state,
		 *            the result will be <code>false</code>.
		 *
		 * @return
		 *            <code>true</code> if the iterator's state changed
		 *            succesfully; otherwise, if the iterator has reached its
		 *            final state, the result will be <code>false</code>.
		 *
		 * @warning
		 *            If not sufficiently memory could be provided, exceptions
		 *            may be thrown correspondingly.
		 */
		bool next();

		/**
		 * @brief
		 *            Access the size of the array from where <i>k</i> entries
		 *            of different positions is selected in each iteration.
		 *
		 * @return
		 *            The size of the array from where <i>k</i> entries
		 *            of different positions is selected in each iteration.
		 */
		inline int get_n() const { return this->n; }

		/**
		 * @brief
		 *            Access the number of entries of different positions that
		 *            are selected from an array of size <i>n</i> in each
		 *            iteration.
		 *
		 * @return
		 *            The number of entries of different positions that are
		 *            selected from an array of size <i>n</i> in each
		 */
		inline int get_k() const { return this->k; }

		/**
		 * @brief
		 *            Access a state of the iterator for internal use.
		 *
		 * @return
		 *            A state of the iterator for internal use.
		 */
		inline int get_idx() const { return this->idx; }

		/**
		 * @brief
		 *            Access the iterator's sub-iterator for internal use.
		 *
		 * @return
		 *			  The iterator's sub-iterator for internal use.
		 */
		inline const BaseBinomialIterator *getSubIteratorPtr() const {
			return this->subIteratorPtr;
		}

	};

	/**
	 * @brief
	 *            Provides a mechanism for iterating through and choosing
	 *            <i>k</i> elements of different positions from arrays of
	 *            size <i>n</i>.
	 *
	 * @details
	 *            Assume we want to iterate through all possible choices
	 *            of selected <i>k</i> integers from an array of size
	 *            <i>n</i>. For example, let <i>n=5</i> and
	 *            <pre>
	 *             int array[] = {5,6,7,8,9};
	 *            </pre>
	 *            Then an iterator that can iterate through all
	 *            choices of <i>k=3</i> entries (of different position) can be
	 *            created via
	 *            <pre>
	 *             BinomialIterator<int> it(5,3);
	 *            </pre>
	 *            Now, if we want to output all choices of <i>k=3</i> entries
	 *            from the array we may run
	 *            <pre>
	 *             do {
	 *             		int selection[3];
	 *             		it.select(selection,array);
	 *
	 *             		for ( int i = 0 ; i < 3 ; i++ ) {
	 *             			cout << selection[i] << " ";
	 *             		}
	 *             		cout << endl;
	 *
	 *             } while ( it.next() );
	 *            </pre>
	 *            whose output will be
	 *            <pre>
	 *            5 6 7
	 *            5 6 8
	 *            5 6 9
	 *            5 7 8
	 *            5 7 9
	 *            5 8 9
	 *            6 7 8
	 *            6 7 9
	 *            6 8 9
	 *            7 8 9
	 *            </pre>
	 *            which corresponds to
	 *            \f[
	 *            \left({n\atop k}\right)=\left({5\atop 3}\right)=10
	 *            \f]
	 *            different choices.
	 */
	template <class T>
	class BinomialIterator : public BaseBinomialIterator {

	public:

		/**
		 * @brief
		 *            Creates an instance of an iterator to iterate through
		 *            all choices of <i>k</i> entries (of different positions)
		 *            from an array of size <i>n</i>.
		 *
		 * @details
		 *            The constructor just calls the constructor of its
		 *            superclass.
		 *
		 * @param n
		 *            Size of the arrays from where <i>k</i> elements are
		 *            chosen in each iteration.
		 *
		 * @param k
		 *            Specifies the number of elements that are chosen from
		 *            an array of size <i>n</i> in each iteration.
		 *
		 * @warning
		 *            If not sufficient memory could be provided, exceptions
		 *            may be thrown correspondingly.
		 *
		 * @warning
		 *            If <i>n</i> or <i>k</i> is smaller than zero or if
		 *            <i>n</i> is smaller than <i>k</i> an error message will
		 *            be written to <code>stderr</code> and an exit with
		 *            status 'EXIT_FAILURE' will be caused.
		 */
		inline BinomialIterator(int n, int k) : BaseBinomialIterator(n,k) {
		}

		/**
		 * @brief
		 *            Selects <i>k</i> elements of different position from
		 *            <code>array</code> corresponding to the iterators state
		 *            and stores the successively in <code>selection</code>.
		 *
		 * @param selection
		 *            Will contain <i>k=</i>\link get_k()\endlink choices
		 *            of different positions from the first
		 *            <i>n</i>=\link get_n()\endlink entries of
		 *            <code>array</code>.
		 *
		 * @param array
		 *            Array from where to choose
		 *            <i>k</i>=\link get_k()\endlink
		 *            entries of different positions according to the
		 *            iterator's current state.
		 *
		 * @warning
		 *            If <code>array</code> does not contain at least <i>n</i>
		 *            valid entries of type <code>T</code> or if
		 *            <code>selection</code> was not allocated successfully to
		 *            hold at least <i>k</i> entries the behavior of the
		 *            method is undocumented. Furthermore, the arrays
		 *            <code>selection</code> and <code>array</code> should not
		 *            cross.
		 */
		void select( T *selection , const T *array ) const;

	};

	/**
	 * @brief
	 *            Selects <i>k</i> elements of different position from
	 *            <code>array</code> corresponding to the iterators state
	 *            and stores the successively in <code>selection</code>.
	 *
	 * @details
	 *            see above.
	 */
	template<class T>
	void BinomialIterator<T>::select( T *selection , const T *array ) const {

		// We follow the iterator and its recursive sub-iterators
		// until we reach 'NULL'.
	    const BaseBinomialIterator *itPtr;

	    // Initialize the pointer by 'this';
	    itPtr = this;

	    // Keeps track of the next index position from where it is possible
	    // to select an array's entry.
	    int k0 = 0;

	    // The next array's entry is set at the 'i'th position of 'selection'.
	    int i = 0;

	    do {

	    	// The 'i'th selection where 'i' is incremented afterwards.
	    	selection[i++] = array[k0+itPtr->get_idx()];

	    	// Next possible index position from where we can choose the
	    	// array's entry.
	    	k0 += itPtr->get_idx()+1;

	    	// Next sub-iterator.
	    	itPtr = itPtr->getSubIteratorPtr();

	    } while ( itPtr != NULL ); // Iteration until the iterator is 'NULL'.
	}

}

#endif /* THIMBLE_BINOMIALITERATOR_H_ */
