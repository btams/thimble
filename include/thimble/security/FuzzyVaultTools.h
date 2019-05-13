/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
 *
 *  Copyright 2013, 2014 Benjamin Tams
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
 * @file FuzzyVaultTools.h
 *
 * @brief
 *            Provides convenience functionalities related with the fuzzy
 *            vault scheme.
 *
 * @details
 *            The fuzzy vault scheme is a cryptographic construction to
 *            protect biometric data.
 *            <br>
 *            This file provides a class that collects utility functions
 *            that are related to the general fuzzy vault scheme. It features
 *            a brute-force attack against the fuzzy vault. Furthermore, a
 *            function that generates a random instance of the fuzzy vault
 *            scheme of specified parameters is provided to allow the user to
 *            perform simulations. A brute-force decoder is provided as well.
 *            For simulations, a function that estimates the brute-force security
 *            of a vault of specified parameters is provided.
 *
 * @author Benjamin Tams
 *
 * @see thimble::FuzzyVaultTools
 */
#ifndef THIMBLE_FUZZYVAULTTOOLS_H_
#define THIMBLE_FUZZYVAULTTOOLS_H_

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
	 *           Provides tools related with the fuzzy vault scheme which is
	 *           to protect biometrics and other noisy data.
	 */
	class THIMBLE_DLL FuzzyVaultTools {

	public:

		/**
		 * @brief
		 *            Creates a random instance of a fuzzy vault with
		 *            specified parameters
		 *
		 * @details
		 *            Using this function, a random fuzzy vault instance
		 *            (defined over the field specified by <code>gf</code>) of
		 *            size <code>n</code> is generated in where a polynomial
		 *            of degree smaller than <code>k</code> interpolates
		 *            exactly <code>t</code> vault points.
		 *            <br><br>
		 *            After invocation, <code>x</code> and <code>y</code>
		 *            contain the vault points <code>(x[i],y[i])</code>
		 *            abscissa and ordinate values, respectively. Furthermore,
		 *            the secret polynomial is returned and can be saved
		 *            if need for simulations
		 *            (e.g., to compute its hash value).
		 *
		 * @param x
		 *            Will contain the vault points abscissa values.
		 *
		 * @param y
		 *            Will contain the vault points ordinate values.
		 *
		 * @param n
		 *            Specifies the size of the vault.
		 *
		 * @param t
		 *            Specifies the number of genuine vault points which
		 *            must be smaller than (or equal) <code>n</code>.
		 *
		 * @param k
		 *            Specifies the size of the secret polynomial
		 *            hidden in the vault which must be smaller than
		 *            (or equal) <code>t</code>.
		 *
		 * @param gf
		 *            The finite field over that the vault is constructed.
		 *
		 * @param tryRandom
		 *            If <code>true</code> the function is advised to use
		 *            a random generator that is cryptographically secure
		 *            (see <code>MathTools::rand8()</code>).
		 *
		 * @return
		 *            The polynomial which is hidden by the created vault.
		 *            The polynomial can be used to check whether a simulated
		 *            attack against the created vault was successful.
		 *
		 * @warning
		 *            In the following case an error message is printed
		 *            to <code>stderr</code> and the program exits with
		 *            status 'EXIT_FAILURE'.
		 *            <ul>
		 *             <li>
		 *              If <code>n</code> is greater than the finite field's
		 *              cardinality.
		 *             </li>
		 *             <li>
		 *              <code>n</code> is smaller than or equal 0.
		 *             </li>
		 *             <li>
		 *              <code>t</code> is smaller than or equal 0.
		 *             </li>
		 *             <li>
		 *              <code>k</code> is smaller than or equal 0.
		 *             </li>
		 *             <li>
		 *              <code>t</code> is greater than <code>n</code>.
		 *             </li>
		 *             <li>
		 *              <code>k</code> is greater than <code>t</code>.
		 *             </li>
		 *             <li>
		 *              Not sufficient memory could be provided to
		 *              create the vault.
		 *             </li>
		 *            </ul>
		 *
		 * @warning
		 *            If <code>x</code> or <code>y</code> cannot store
		 *            at least <code>n</code> unsigned 32-bit integers,
		 *            the function may run into undocumented behavior.
		 *
		 */
		static SmallBinaryFieldPolynomial createRandomInstance
		( uint32_t *x , uint32_t *y , int n , int t , int k ,
		  const SmallBinaryField & gf , bool tryRandom = true );

        /**
         * @brief
         *            Creates a random instance of the improved fuzzy vault scheme.
         *
         * @details
         *            The function generates a random polynomial <i>f</i> of degree
         *            smaller <i>k</i>; furthermore, a product
         *            \f[
         *             \prod_{j=0}^{t-1}(X-x_j)
         *            \f]
         *            where the \f$x_j\f$ are selected to be random and distinct.
         *            Then the function outputs
         *            \f[
         *             V(X)=f(X)+\prod_{j=0}^{t-1}(X-x_j).
         *            \f]
         *            The function also returns <i>f</i> to enable evaluations of
         *            attacks in order to decide whether they were successful or
         *            not.
         *
         * @param V
         *            The random vault polynomial.
         *
         * @param t
         *            The number of random genuine features \f$x_j\f$ protected by
         *            the vault.
         *
         * @param k
         *            The length of the vault's secret polynomial.
         *
         * @param tryRandom
         *            If <code>true</code>, the function is advised to use a
         *            cryptographic number generator; otherwise, the method
         *            wraps around the standard <code>rand()</code> function.
         *
         * @warning
         *            If the relation \f$0\leq k<t\leq n\f$, where \f$n\f$ is
         *            the cardinality of the field over which \f$V\f$ is
         *            defined, the function prints an error message to
         *            <code>stderr</code> and exits with status
         *            'EXIT_FAILURE'.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        static SmallBinaryFieldPolynomial createRandomInstance
        ( SmallBinaryFieldPolynomial & V , int t , int k ,
          bool tryRandom = true );

		/**
		 * @brief
		 *            Attempts to break an instance of the fuzzy vault scheme.
		 *
		 * @details
		 *            A vault of size <code>n</code> is passed to the method
		 *            where its points are <code>(x[i],y[i])</code> for
		 *            <code>i=0,...,n-1</code>. The field over that the vault
		 *            is defined is the same as the field over that the
		 *            polynomial <code>f</code> was initialized.
		 *            <br><br>
		 *            In each iteration, the attack selects <code>k</code>
		 *            vault points at random (using the standard
		 *            <code>rand()</code> function), determines its
		 *            interpolation polynomial, and checks whether its 160-bit
		 *            SHA-1 hash value agrees with the value contained in
		 *            <code>hash</code>. If successful, <code>f</code> will be
		 *            assigned with the found polynomial and the function
		 *            returns <code>true</code>. Otherwise, the iteration
		 *            continues until a polynomial with the given hash is
		 *            found or if <code>maxIts</code> iterations have been
		 *            performed. If no polynomial could be found this way,
		 *            the function returns <code>false</code> and leaves
		 *            <code>f</code> unchanged.
		 *
		 * @param f
		 *            Will contain the secret vault polynomial if the function
		 *            returns <code>true</code>.
		 *
		 * @param x
		 *            Contains the <code>n</code> successive vault point's
		 *            abscissas.
		 *
		 * @param y
		 *            Contains the <code>n</code> successive vault point's
		 *            ordinate values.
		 *
		 * @param n
		 *            The vault's size.
		 *
		 * @param k
		 *            The size of the secret polynomial (i.e. the polynomial
		 *            is assumed to be of degree smaller than <code>k</code>)
		 *
		 * @param hash
		 *            The SHA-1 hash value of the secret polynomial.
		 *
		 * @param maxIts
		 *            The maximal number of iterations that are performed
		 *            before the attack stops.
		 *
		 * @warning
		 *            In the following cases the function prints an error
		 *            message to <code>stderr</code> and exits with status
		 *            -1.
		 *            <ul>
		 *             <li>
		 *              <code>n</code> is smaller than or equal 0
		 *              or if <code>n</code> is larger than RAND_MAX.
		 *             </li>
		 *             <li>
		 *              <code>k</code> is smaller than or equal 0
		 *              or if <code>k</code> is larger than <code>n</code>.
		 *             </li>
		 *            </ul>
		 *
		 * @warning
		 *            In the following cases the function may run into
		 *            unexpected, undocumented behavior.
		 *            <ul>
		 *             <li>
		 *              If <code>x</code> or <code>y</code> do not contain
		 *              at least <code>n</code> valid elements in the finite
		 *              field <code>f.getField()</code>.
		 *             </li>
		 *             <li>
		 *              If <code>hash</code> does not contain 5 valid
		 *              unsigned 32-bit integers.
		 *             </li>
		 *            </ul>
		 */
		static bool bfattack
		( SmallBinaryFieldPolynomial & f ,
		  const uint32_t *x , const uint32_t *y ,
		  int n , int k , const uint32_t hash[5] ,
		  uint64_t maxIts );

		/**
		 * @brief
		 *            Attempts to break an instance of the fuzzy vault scheme.
		 *
		 * @details
		 *            A vault of size <code>n</code> is passed to the method
		 *            where its points are <code>(x[i],y[i])</code> for
		 *            <code>i=0,...,n-1</code>. The field over that the vault
		 *            is defined is the same as the field over that the
		 *            polynomial <code>f</code> was initialized.
		 *            <br><br>
		 *            In each iteration, the attack selects <code>k</code>
		 *            vault points at random (using the standard
		 *            <code>rand()</code> function), determines its
		 *            interpolation polynomial, and checks whether its 160-bit
		 *            SHA-1 hash value agrees with the value contained in
		 *            <code>hash</code>. If successful, <code>f</code> will be
		 *            assigned with the found polynomial and the function
		 *            returns <code>true</code>. Otherwise, the iteration
		 *            continues until a polynomial with the given hash is
		 *            found or if <code>maxIts</code> iterations have been
		 *            performed. If no polynomial could be found this way,
		 *            the function returns <code>false</code> and leaves
		 *            <code>f</code> unchanged.
		 *
		 * @param f
		 *            Will contain the secret vault polynomial if the function
		 *            returns <code>true</code>.
		 *
		 * @param x
		 *            Contains the <code>n</code> successive vault point's
		 *            abscissas.
		 *
		 * @param y
		 *            Contains the <code>n</code> successive vault point's
		 *            ordinate values.
		 *
		 * @param n
		 *            The vault's size.
		 *
		 * @param k
		 *            The size of the secret polynomial (i.e. the polynomial
		 *            is assumed to be of degree smaller than <code>k</code>)
		 *
		 * @param hash
		 *            The SHA-1 hash value of the secret polynomial.
		 *
		 * @param maxIts
		 *            The maximal number of iterations that are performed
		 *            before the attack stops.
		 *
		 * @warning
		 *            In the following cases the function prints an error
		 *            message to <code>stderr</code> and exits with status
		 *            -1.
		 *            <ul>
		 *             <li>
		 *              <code>n</code> is smaller than or equal 0
		 *              or if <code>n</code> is larger than RAND_MAX.
		 *             </li>
		 *             <li>
		 *              <code>k</code> is smaller than or equal 0
		 *              or if <code>k</code> is larger than <code>n</code>.
		 *             </li>
		 *            </ul>
		 *
		 * @warning
		 *            In the following cases the function may run into
		 *            unexpected, undocumented behavior.
		 *            <ul>
		 *             <li>
		 *              If <code>x</code> or <code>y</code> do not contain
		 *              at least <code>n</code> valid elements in the finite
		 *              field <code>f.getField()</code>.
		 *             </li>
		 *             <li>
		 *              If <code>hash</code> does not contain 20 valid
		 *              unsigned 8-bit integers.
		 *             </li>
		 *            </ul>
		 */
		static bool bfattack
		( SmallBinaryFieldPolynomial & f ,
		  const uint32_t *x , const uint32_t *y ,
		  int n , int k , const uint8_t hash[20] ,
		  uint64_t maxIts );

		/**
		 * @brief
		 *           Attempts to decode a polynomial given unlocking points.
		 *
		 * @details
		 *           Let \f${\bf U}=\{(x[0],y[0]),...,(x[t-1],y[t-1])\}\f$ be
		 *           the unlocking set. Then the decoder iterates through all
		 *           possible choices of <i>k</i> subsets of the unlocking
		 *           sets. Then its corresponding candidate polynomial
		 *           \f$f^*(X)\f$ of degree smaller <i>k</i> is computed which
		 *           interpolates the points in the subset. The SHA-1 hash value
		 *           of the candidate polynomial is then compared with the
		 *           hash value given by <code>hash</code>; if both are equal,
		 *           the function assigns <i>f</i> with the candidate polynomial
		 *           and returns <code>true</code>; otherwise, if no candidate
		 *           polynomial with the hash of the correct polynomial can be
		 *           found, <i>f</i> will be left unchanged and the function
		 *           returns <code>false</code>.
		 *
		 * @param f
		 *           If decoding was successful, the polynomial will be equals
		 *           the decoded polynomial.
		 *
		 * @param x
		 *           Successive abscissas of the unlocking set.
		 *
		 * @param y
		 *           Successive ordinates of the unlocking set.
		 *
		 * @param n
		 *           Size of the unlocking set.
		 *
		 * @param k
		 *           Size of the secret polynomial.
		 *
		 * @param hash
		 *           The SHA-1 hash value of the secret polynomial.
		 * @return
		 *           <code>true</code> if decoding was successful and
		 *           <code>false</code> otherwise.
		 *
		 * @warning
		 *           If the first <i>t</i> entries of <i>x</i>
		 *           are not pairwise distinct, an error message will be
		 *           printed to <code>stderr</code> and causes an exit
		 *           with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *           If <i>x</i> or if <i>y</i> do not contain at least
		 *           <i>t</i> valid elements in the finite field over which
		 *           <i>f</i> is defined, the function runs into unexpected
		 *           behavior.
		 *
		 * @warning
		 *           If <i>n</i> or <i>k</i> are smaller than zero
		 *           the function prints an error message to
		 *           <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *           If not sufficient memory can be allocated to run the
		 *           function the function prints an error message to
		 *           <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		static bool bfdecode
		( SmallBinaryFieldPolynomial & f ,
		  const uint32_t *x , const uint32_t *y ,
		  int n , int k , const uint32_t hash[5] );

		/**
		 * @brief
		 *           For a random fuzzy vault of specified parameters,
		 *           the function computes the difficulty for choosing
		 *           the <i>k</i> genuine points in the vault.
		 *
		 * @details
		 *           More precisely, the function computes
		 *           \f[
		 *            {\bf bf}(n,t,k):=\left({n\atop k}\right)\cdot
		 *            \left({t\atop k}\right)^{-1}
		 *           \f]
		 *           such that \f${\bf bf}(n,t,k)^{-1}\f$ is the probability
		 *           that <i>k</i> randomly selected vault points will be
		 *           all genuine points.
		 *
		 * @param n
		 *           Size of the vault.
		 *
		 * @param t
		 *           Number of genuine points in the vault.
		 *
		 * @param k
		 *           Size of secret polynomial hidden in the vault.
		 *
		 * @return
		 *           The difficulty in choosing <i>k</i> genuine points
		 *           from the vault points.
		 */
		static long double bfsecurity( int n , int t , int k );
	};
}


#endif /* THIMBLE_FUZZYVAULTTOOLS_H_ */
