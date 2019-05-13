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
 * @file EEAAttack.h
 *
 * @brief
 *            Provides a class for running a record multiplicity attack
 *            against the <em>improved fuzzy vault scheme</em>.
 *
 * @details
 * @section sec_arm_improved_fv Attack via Record Multiplicity against the Improved Fuzzy Vault Scheme
 *
 * In this tutorial, we show how multiple instances of the improved
 * fuzzy vault scheme can be cross-matched and, if successful, even
 * (in parts and sometimes even fully) be broken. The theory of the
 * attack can be consulted in the following reference
 * <ul>
 *  <li>
 *   <b>Merkle and Tams (2013).</b> Security of the Improved Fuzzy Vault
 *   Scheme in the Presence of Record Multiplicity ,
 *   <a href="http://arxiv.org/abs/1312.5225" target="_blank">arXiv:1312.5225</a>.
 *  </li>
 * </ul>
 * and in this tutorial, it is shown how the \link thimble::EEAAttack
 * EEAAttack\endlink class from THIMBLE can be used to reproduce
 * the experimental result given in the above reference.
 *
 * An instance of the improved fuzzy vault scheme defined over the
 * field
 * <pre>
 *  SmallBinaryField gf(16);
 * </pre>
 * say, can be created with the help of
 * the \link thimble::SmallBinaryFieldPolynomial
 * SmallBinaryFieldPolynomial\endlink class. Therefore, the parameters
 * of the improved fuzzy vault scheme should be selected first in which
 * feature sets of size
 * <pre>
 *    int t = 44;
 * </pre>
 * (say) are protected using a secret polynomial of size
 * <pre>
 *    int k = 10;
 * </pre>
 * (say). Now, let
 * <pre>
 *  uint32_t *A , B*;
 * </pre>
 * be two arrays of length <code>t</code> that contain exactly
 * <pre>
 *    int omega = 30;
 * </pre>
 * (say) common integers of type <code>uint32_t</code>
 * representing elements from the galois field <code>gf</code>.
 *
 * Each of the feature sets <code>A</code> and <code>B</code> are binded
 * to polynomials of size <code>k</code> being selected randomly using
 * \link thimble::SmallBinaryFieldPolynomial::random()
 * SmallBinaryFieldPolynomial::random()\endlink:
 * <pre>
 *    SmallBinaryFieldPolynomial f(gf) , g(gf);
 *
 *    f.random(k,true);
 *    g.random(k,true);
 * </pre>
 * To build the vaults, the characteristic polynomials, e.g., computed
 * with the help
 * of \link thimble::SmallBinaryFieldPolynomial::buildFromRoots()
 * SmallBinaryFieldPolynomial::buildFromRoots()\endlink
 * <pre>
 *    SmallBinaryFieldPolynomial chiA(gf), chiB(gf);
 *
 *    chiA.buildFromRoots(A,t);
 *    chiB.buildFromRoots(B,t);
 * </pre>
 * are added to the secret polynomials:
 * <pre>
 *    SmallBinaryFieldPolynomial V(gf) , W(gf);
 *
 *    V = f + chiA;
 *    W = g + chiB;
 * </pre>
 *
 * @subsection sec_eeaattack_crossmatching Cross-Matching
 *
 * In the following, we assume that we are an attacker who
 * only has access to the vaults <code>V</code> and <code>W</code>
 * but knows nothing about <code>A</code> or <code>B</code>.
 * Using the \link thimble::EEAAttack EEAAttack\endlink we can
 * try to recover the differences between the feature sets
 * <code>A</code> and <code>B</code> via
 * <pre>
 *    EEAAttack attack;
 *
 *    bool success = attack.perform(V,W,k);
 *
 *    if ( success ) {
 *       cout << "RELATED" << endl;
 *    } else {
 *       cout << "NON-RELATED" << endl;
 *    }
 * </pre>
 * which outputs
 * <pre>
 *    RELATED
 * </pre>
 * whenever the relation
 * <pre>
 *    omega >= (t+k)/2;
 * </pre>
 * is fulfilled. If <code>omega<(t+k)/2</code>, then the attack can
 * be successful and falsely label two vault as related. However,
 * the probability that this happens may be very small for parameters that
 * we expect to encounter in practice. Thus, even though the relation between
 * two instances of the improved fuzzy vault scheme can be examined quite
 * efficiently with the attack, i.e., instances of the improved fuzzy vault
 * scheme are vulnerable to cross-matching
 *
 * @subsection sec_eeaattack_reversibility Reversibility
 *
 * If \link thimble::EEAAttack::perform() attack.perform(V,W,k)\endlink
 * returns <code>true</code> and if the difference between the feature sets
 * <code>A</code> and <code>B</code> is larger than <code>k</code>, then
 * the vaults can even be fully broken. Therefore, we may run
 * <pre>
 *    SmallBinaryFieldPolynomial q(gf);
 *    q.buildFromRoots(attack.getFeatures1(),attack.getNumFeatures1());
 *
 *    SmallBinaryFieldPolynomial rec_f = V % q;
 * </pre>
 * where <code>rec_f</code> equals the secret polynomial <code>f</code>
 * if <code>attack.getNumFeatures1()>=k</code>. Analogously, the polynomial
 * <code>g</code> can be recovered via
 * <pre>
 *    SmallBinaryFieldPolynomial p(gf);
 *    p.buildFromRoots(attack.getFeatures2(),attack.getNumFeatures2());
 *
 *    SmallBinaryFieldPolynomial rec_g = W % p;
 * </pre>
 * if <code>attack.getNumFeatures2()>=k</code>. After <code>f</code>
 * or <code>g</code> have been recovered successfully, the feature sets
 * <code>A</code> and <code>B</code> can even be recovered from the vaults
 * by finding the roots of the polynomial <code>V-rec_f</code> and
 * <code>W-rec_g</code>, respectively, using
 * the \link thimble::SmallBinaryFieldPolynomial::findRoots()
 * SmallBinaryFieldPolynomial::findRoots()\endlink function.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_EEAATTACK_H_
#define THIMBLE_EEAATTACK_H_

#include <stdint.h>
#include <vector>

#include <thimble/dllcompat.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Instances of this class provide a function for
	 *            running an extended Euclidean algorithm-based record
	 *            multiplicity attack against the improved fuzzy
	 *            vault scheme and store the result.
	 */
	class THIMBLE_DLL EEAAttack {

	public:

		/**
		 * @brief
		 *            Standard constructor.
		 */
		EEAAttack();

		/**
		 * @brief
		 *            Destructor.
		 *
		 * @details
		 *            Frees all memory that has been allocated for
		 *            the members of this class.
		 */
		~EEAAttack();

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
		bool perform
		( const SmallBinaryFieldPolynomial & V ,
		  const SmallBinaryFieldPolynomial & W , int k );

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
		int getNumOverlap() const;

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
		int getNumFeatures1() const;

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
		int getNumFeatures2() const;

		/**
		 * @brief
		 *            An array that forms a candidate set for the explicit
		 *            elements that are protected by the first vault
		 *            after the attack \link perform()\endlink has been
		 *            run successfully and that are not protected by
		 *            the second vault.
		 *
		 * @return
		 *            An array containing <code>getNumFeatures1()</code>
		 *            integers of type <code>uint32_t</code>.
		 */
		const uint32_t *getFeatures1() const;

		/**
		 * @brief
		 *            An array that forms a candidate set for the explicit
		 *            elements that are protected by the second vault
		 *            after the attack \link perform()\endlink has been
		 *            run successfully and that are not protected by
		 *            the first vault.
		 *
		 * @return
		 *            An array containing <code>getNumFeatures2()</code>
		 *            integers of type <code>uint32_t</code>.
		 */
		const uint32_t *getFeatures2() const;

		/**
		 * @brief
		 *            Clears all data stored by this object such that
		 *            its state is equivalent as after construction.
		 *
		 * @see EEAAttack()
		 */
		void clear();

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
		void swap();

	private:

		/**
		 * @brief
		 *            A candidate for the number of common feature
		 *            set elements protected by two vaults for which the
		 *            attack has been applied.
		 *
		 * @see getNumOverlap()
		 */
		int num_overlap;

		/**
		 * @brief
		 *            The number of feature elements stored
		 *            in \link features1\endlink.
		 *
		 * @see getNumFeatures1()
		 */
		int num_features1;

		/**
		 * @brief
		 *            The number of feature elements stored
		 *            in \link features2\endlink.
		 *
		 * @see getNumFeatures2()
		 */
		int num_features2;

		/**
		 * @brief
		 *            Contains candidates of the feature elements protected
		 *            by the first vault for which this attack has been
		 *            applied that are not protected by the second vault.
		 *
		 * @details
		 *            This array holds \link num_features1\endlink distinct
		 *            elements of type <code>uint32_t</code>.
		 *
		 * @see getFeatures1()
		 */
		uint32_t *features1;

		/**
		 * @brief
		 *            Contains candidates of the feature elements protected
		 *            by the second vault for which this attack has been
		 *            applied that are not protected by the first vault.
		 *
		 * @details
		 *            This array holds \link num_features2\endlink distinct
		 *            elements of type <code>uint32_t</code>.
		 *
		 * @see getFeatures2()
		 */
		uint32_t *features2;
	};
}

#endif /* THIMBLE_EEAATTACK_H_ */
