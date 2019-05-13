/*
 *  THIMBLE --- Research Library for Development and Analysis of
 *  Fingerprint-Based Biometric Cryptosystems.
 *
 *  Copyright 2014 Benjamin Tams
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
 * @file EEAAttack.cpp
 *
 * @brief
 *            Implementation of a class that runs a record multiplicity
 *            attack against the <em>improved fuzzy vault scheme</em>
 *            as provided by the 'EEAAttack.h' header.
 *
 * @author Benjamin Tams
 */

#include <cmath>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>
#include <thimble/security/EEAAttack.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Standard constructor.
	 */
	EEAAttack::EEAAttack() {
		this->num_overlap         = -1;
		this->num_features1       = 0;
		this->num_features2       = 0;
		this->features1           = NULL;
		this->features2           = NULL;
	}

	/**
	 * @brief
	 *            Destructor.
	 *
	 * @details
	 *            Frees all memory that has been allocated for
	 *            the members of this class.
	 */
	EEAAttack::~EEAAttack() {

		this->clear();
	}

	/**
	 * @brief
	 *            Given two instances of the improved fuzzy vault
	 *            scheme, this function attempts to determine whether
	 *            they are related and, if labeled as related, recover
	 *            the differences between the protected features
	 *            explicitly.
	 *
	 * @param V
	 *            The first instance of the improved fuzzy vault
	 *            scheme.
	 *
	 * @param W
	 *            The second instance of the improved fuzzy vault
	 *            scheme.
	 *
	 * @param k
	 *            The (maximal) size of the secret polynomials protected
	 *            by the two vaults <code>V</code> and <code>W</code>.
	 *
	 * @return
	 *            <code>true</code> if the two vaults were positively
	 *            labeled as being related in which case candidates
	 *            for the protected feature sets differences are stored
	 *            by this object; otherwise, this function returns
	 *            <code>false</code> in which any previously computed
	 *            data stored by this object is cleared.
	 *
	 * @warning
	 *           If not enough memory could be provided, an error
	 *           message is printed to <code>stderr</code> and the
	 *           program exits with status 'EXIT_FAILURE'.
	 */
	bool EEAAttack::perform
	( const SmallBinaryFieldPolynomial & V ,
	  const SmallBinaryFieldPolynomial & W , int k ) {

		int t , s;
		t = V.deg();
		s = W.deg();

		if ( t < s ) {

			// In the attack, we assume that the first vault protects
			// the larger feature set. If this assumption is violated,
			// we run the attack with both vaults exchanged ...
			bool state = perform(W,V,k);
			// ... which requires to swap the candidates for the
			// recovered feature sets' differences (if any).
			this->swap();

			return state;
		}

		// Clear any previously result from the attack.
		this->clear();

		// Run the extended Euclidean algorithm ...
		vector<SmallBinaryFieldPolynomial> R , P , Q;
		xgcd(R,P,Q,V,W);

		// .. and determine the relation that minimizes
		// epsilon = deg(Q[j0]) where ...
		int j0 = -1 , minEpsilon = INT_MAX;
		for ( int j = 0 ; j < (int)R.size() ; j++ ) {
			// ... R[j0] != 0  and deg(Q[j0])+k > R[j0].
			if ( !R[j].isZero() && R[j].deg() < k + Q[j].deg() ) {
				if ( Q[j].deg() < minEpsilon ) {
					minEpsilon = Q[j].deg();
					j0 = j;
				}
			}
		}

		// If none such j0 exist, the two vaults definitely do not
		// protect feature sets that overlap in at least (t+k)/2
		// elements and they are labeled as non-related.
		if ( j0 < 0 ) {
			return false;
		}

		// Runs Step 3) of the algorithm.
		SmallBinaryFieldPolynomial tmp(V.getField());
		rem(tmp,V,Q[j0]);
		if ( tmp.deg() >= k ) {
			return false;
		}

		// Now, the roots of Q[j0] and P[j0] form the
		// differences of A\B and B\A, respectively, where
		// A denotes the feature set protected by V and
		// B the feature set protected by W.

		// Allocate memory such that 'feature1' can store
		// the elements of A\B.
		if ( Q[j0].deg() == 0 ) {
			this->features1 = NULL;
		} else {
			this->features1 = (uint32_t*)malloc( Q[j0].deg() * sizeof(uint32_t) );
			if ( this->features1 == NULL ) {
				cerr << "PartialRecoveryAttack::perform: out of memory." << endl;
				exit(EXIT_FAILURE);
			}
		}

		// Allocate memory such that 'feature2' can store
		// the elements of B\A.
		if ( P[j0].deg() == 0 ) {
			this->features2 = NULL;
		} else {
			this->features2 = (uint32_t*)malloc( P[j0].deg() * sizeof(uint32_t) );
			if ( this->features2 == NULL ) {
				cerr << "PartialRecoveryAttack::perform: out of memory." << endl;
				exit(EXIT_FAILURE);
			}
		}

		// Find the roots of P[j0].
		this->num_features2 = P[j0].findRoots(this->features2);
		if ( this->num_features2 != P[j0].deg() ) {
			// If P[j0] does not completely split into linear factors
			// the feature set protected by W that overlaps the feature set
			// protected by V does not exist in the base field (but
			// possibly in an extension). In any case, the two vaults
			// are then labled as non-related.
			clear();
			return false;
		}

		// Find the roots of Q[j0].
		this->num_features1 = Q[j0].findRoots(this->features1);
		if ( this->num_features1 != Q[j0].deg() ) {
			// see above; the same as for P[j0]
			clear();
			return false;
		}

		// Finally, determine the number of common feature set
		// elements which might exist in an extension field only of which
		// probability, depending on the parameters, can be very small.
		this->num_overlap = t - this->num_features1;

		// Label the two vaults as related.
		return true;
	}

	/**
	 * @brief
	 *            Returns a candidate of the number of common feature
	 *            elements protected by two vault against which the
	 *            attack \link perform()\endlink has been successfully
	 *            applied.
	 *
	 * @return
	 *            Candidate for the number of common feature elements
	 *            if the attack \link perform()\endlink was previously
	 *            applied successfully; otherwise, the result of this
	 *            function is -1.
	 */
	int EEAAttack::getNumOverlap( ) const {
		return this->num_overlap;
	}

	/**
	 * @brief
	 *            Returns a candidate for the number of feature elements
	 *            protected by the first vault against which the
	 *            attack \link perform()\endlink has been successfully
	 *            applied and that are not protected by the second
	 *            vault.
	 *
	 * @return
	 *            If the attack \link perform()\endlink returned
	 *            <code>true</code>, this function returns the candidate
	 *            for the number that are protected by the first vault
	 *            and not by the second; otherwise,
	 *            if \link perform()\endlink was <code>false</code>, this
	 *            function returns <code>0</code>.
	 */
	int EEAAttack::getNumFeatures1() const {
		return this->num_features1;
	}

	/**
	 * @brief
	 *            Returns a candidate for the number of feature elements
	 *            protected by the second vault against which the
	 *            attack \link perform()\endlink has been successfully
	 *            applied and that are not protected by the first
	 *            vault.
	 *
	 * @return
	 *            If the attack \link perform()\endlink returned
	 *            <code>true</code>, this function returns the candidate
	 *            for the number that are protected by the second vault
	 *            and not by the first; otherwise,
	 *            if \link perform()\endlink was <code>false</code>, this
	 *            function returns <code>0</code>.
	 */
	int EEAAttack::getNumFeatures2() const {
		return this->num_features2;
	}

	const uint32_t *EEAAttack::getFeatures1() const {
		return this->features1;
	}

	const uint32_t *EEAAttack::getFeatures2( ) const {
		return this->features2;
	}

	/**
	 * @brief
	 *            Clears all data stored by this object such that
	 *            its state is equivalent as after construction.
	 *
	 * @see EEAAttack()
	 */
	void EEAAttack::clear() {
		this->num_overlap = -1;
		this->num_features1 = 0;
		this->num_features2 = 0;
		free(this->features1);
		free(this->features2);
		this->features1 = NULL;
		this->features2 = NULL;
	}

	/**
	 * @brief
	 *            Exchanges the two candidates for the feature
	 *            sets' differences (if non-empty).
	 *
	 * @details
	 *            After swapped, the result
	 *            of \link getNumFeatures1()\endlink
	 *            and \link getNumFeatures2()\endlink will be exchanged
	 *            as well as the result
	 *            of \link getFeatures1()\endlink
	 *            and \link getFeatures2()\endlink
	 */
	void EEAAttack::swap() {

		uint32_t *features;
		int num_features;

		features            = this->features1;
		num_features        = this->num_features1;

		this->features1     = this->features2;
		this->num_features1 = this->num_features2;

		this->features2     = features;
		this->num_features2 = num_features;
	}

}



