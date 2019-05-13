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
 * @file BigInteger.cpp
 *
 * @brief
 *            Implements the functionalities that are provided by the
 *            'BigInteger.h' header which Provides functionalities for
 *            arithmetic of integers of arbitrary precision.
 *
 * @author Benjamin Tams
 *
 * @see thimble::BigInteger
 */
#include <stdint.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include <thimble/math/MathTools.h>
#include <thimble/math/numbertheory/BigInteger.h>

using namespace std;

/**
 * @brief The library's namespace
 */
namespace thimble {

	/**
	 * @brief
	 *            Low-level function comparing two positive big integers.
	 *
	 * @details
	 *            The first big integer is given by the first <i>m</i>
	 *            elements in <i>a</i> which is
	 *            \[
	 *             \sum_{i=0}^{m-1}a[i]\cdot 2^{i\cdot 32}.
	 *            \]
	 *            Similarly, the second big integer is
	 *            \[
	 *             \sum_{i=0}^{n-1}b[i]\cdot 2^{i\cdot 32}.
	 *            \]
	 *            The function returns -1 if the first integer is smaller than
	 *            the second; it returns 1 if the first integer is greater than
	 *            the second; otherwise, if both are equal, it returns 0.
	 *
	 * @param a
	 *            Contains the data for the first big integer.
	 *
	 * @param m
	 *            Specifies the number of relevant <code>uint32_t</code>s
	 *            that define the first big integer.
	 *
	 * @param a
	 *            Contains the data for the second big integer.
	 *
	 * @param m
	 *            Specifies the number of relevant <code>uint32_t</code>s
	 *            that define the second big integer.
	 *
	 * @return
	 *            -1 if the first big integer is smaller than the second; 1 if
	 *            the first big integer is greater than the second; otherwise,
	 *            if both are equal, the function returns 0.
	 *
	 * @warning
	 *            If <i>m</i> or <i>n</i> are smaller than or equals 0 or if
	 *            <i>a</i> or <i>b</i> do not contain at leas <i>m</i> or
	 *            <i>n</i> valid <code<uint32_t</code>, respectively, the
	 *            behavior of the function is undocumented.
	 */
	static int compare
	( register uint32_t *a , register int m ,
	  register uint32_t *b , register int n ) {


		// Ensure 'a' and 'm' correspond to the integer with more elements
		// than 'b' and 'n'. If correspondences have to be swapped therefore,
		// 's' will be set to -1 and 1 otherwise.
		int s;
		if ( m < n ) {

			uint32_t *c;
			int l;

			c = b;
			l = n;

			b = a;
			m = n;

			a = c;
			m = l;

			s = -1;

		} else {
			s = 1;
		}

		// If any of the coefficients in 'a' that are of larger index than
		// 'n' is non-zero the integer defined by 'a' and 'm' must be larger.
		// Because correspondences could have been swapped before, the result
		// is 's' which would be -1 on swapping.
		register int i;
		for ( i = m-1 ; i >= n ; i-- ) {
			if ( a[i] ) {
				return s;
			}
		}

		// Here, 'm' and 'n' are equal: Then the largest non-zero leading
		// coefficient of the determines which one is the largest/smallest
		// integer. Again, if the correspondences were swapped before,
		// the result is depends on 's' also.
		for ( i = n-1 ; i >= 0 ; i-- ) {
			if ( a[i] > b[i] ) {
				return s;
			} else if ( a[i] < b[i] ) {
				return -s;
			}
		}

		// No larger or smaller integer: Both must be equal.
		return 0;
	}

	/**
	 * @brief
	 *             Low-level function that forms the sum of two nonnegative
	 *             multi-precision integers.
	 *
	 * @details
	 *             The integer
	 *             \f[
	 *              \sum_{i=0}^{m-1}a[i]\cdot 2^{i\cdot 32}
	 *             \f]
	 *             is added to the integer
	 *             \f[
	 *              \sum_{i=0}^{n-1}b[i]\cdot 2^{i\cdot 32}
	 *             \f]
	 *             and stored in the array <i>c</i> such that
	 *             \f[
	 *         carry\cdot2^{l\cdot 32}+\sum_{i=0}^{l-1}c[i]\cdot 2^{i\cdot 32}
	 *             \f]
	 *             is the sum. Thereby \f$l=\max(m,n)\f$ and <i>carry</i> is
	 *             the returned result of the function.
	 *             <br><br>
	 *             <b>Note:</b> The arrays <i>c</i>, <i>a</i>, and <i>b</i>
	 *             are allowed to be pairwise equal.
	 *
	 * @param c
	 *             Will contain the coefficients of the sum.
	 *
	 * @param a
	 *             Contains the coefficients of the first integer.
	 *
	 * @param m
	 *             Specifies the number of coefficients for the first integer.
	 *
	 * @param b
	 *             Contains the coefficients of the second integer.
	 *
	 * @param n
	 *             Specifies the number of coefficients of the second integer.
	 *
	 * @return
	 *             The carry of the sum as described in the details.
	 *
	 * @warning
	 *             In the following cases the behavior of the function is
	 *             undocumented.
	 *             <ul>
	 *              <li>
	 *               <i>a</i> does not contain at least <i>m</i> valid
	 *               <code>uint32_t</code>s.
	 *              </li>
	 *              <li>
	 *               <i>b</i> does not contain at least <i>m</i> valid
	 *               <code>uint32_t</code>s.
	 *              </li>
	 *              <li>
	 *               <i>c</i> was not allocated to hold at least
	 *               <i>l=max(m,n)</i> valid <code>uint32_t</code>s.
	 *              </li>
	 *              <li>
	 *               <i>m</i> or <i>n</i> is negative.
	 *              </li>
	 *              <li>
	 *               If the specified arrays overlap when they are not pairwise
	 *               equal.
	 *              </li>
	 *             </ul>
	 */
	static uint32_t add
	( register uint32_t *c ,
	  register const uint32_t *a , register int m ,
	  register const uint32_t *b , register int n ) {

		register uint32_t carry = 0;
		register uint64_t e;
		register int l;
		register int i;

		// How many coefficients to add up?
		l = m<n?m:n;

		// Add the first min(m,n) coefficients
		for ( i = 0 ; i < l ; i++ , a++ , b++ , c++ ) {
			e = (uint64_t)(*a)+(uint64_t)(*b)+(uint64_t)carry;
			*c = (e&0xFFFFFFFF);
			carry = (e>>32);
		}

		// The loop will be entered only if 'm>n': In this
		// case, we add the carry to the elements in 'a'
		// and store the result in 'c'.
		for ( ; i < m ; i++ , a++ , c++ ) {
			e = (uint64_t)(*a) + (uint64_t)carry;
			*c = (e&0xFFFFFFFF);
			carry = (e>>32);
		}

		// The loop will be entered only if 'm<n': In this
		// case, we add the carry to the elements in 'b'
		// and store the result in 'c'.
		for ( ; i < n ; i++ , b++ , c++ ) {
			e = (uint64_t)(*b) + (uint64_t)carry;
			(*c) = (e&0xFFFFFFFF);
			carry = (e>>32);
		}

		return carry;
	}

	/**
	 * @brief
	 *             Low-level function that Forms the difference between two
	 *             nonnegative multi-precision integers.
	 *
	 * @details
	*              The integer
	 *             \f[
	 *              \sum_{i=0}^{n-1}b[i]\cdot 2^{i\cdot 32}
	 *             \f]
	 *             is subtracted from the integer
	 *             \f[
	 *              \sum_{i=0}^{m-1}a[i]\cdot 2^{i\cdot 32}
	 *             \f]
	 *             and stored in the array <i>c</i> such that
	 *             \f[
	 *              2^{l\cdot 32}+\sum_{i=0}^{n-1}c[i]\cdot 2^{i\cdot 32}
	 *             \f]
	 *             is the difference.
	 *             <br><br>
	 *             It is required that the length of <i>a</i> is larger than
	 *             or equals the length of <i>b</i>, i.e. <i>m>=n</i>.
	 *             Moreover, the minuend must be larger than or equals the
	 *             subtrahend.
	 *             <br><br>
	 *             <b>Note:</b> The arrays <i>c</i>, <i>a</i>, and <i>b</i>
	 *             are allowed to be pairwise equal.
	 *
	 * @param c
	 *             Will contain the coefficients of the difference
	 *
	 * @param a
	 *             Contains the coefficients of the minuend.
	 *
	 * @param m
	 *             Specifies the number of coefficients for the minuend.
	 *
	 * @param b
	 *             Contains the coefficients of the subtrahend.
	 *
	 * @param n
	 *             Specifies the number of coefficients of subtrahend.
	 *
	 * @warning
	 *             In the following cases the behavior of the function is
	 *             undocumented.
	 *             <ul>
	 *              <li>
	 *               <i>a</i> does not contain at least <i>m</i> valid
	 *               <code>uint32_t</code>s.
	 *              </li>
	 *              <li>
	 *               <i>b</i> does not contain at least <i>m</i> valid
	 *               <code>uint32_t</code>s.
	 *              </li>
	 *              <li>
	 *               <i>c</i> was not allocated to hold at least
	 *               <i>m</i> valid <code>uint32_t</code>s.
	 *              </li>
	 *              <li>
	 *               <i>m</i> or <i>n</i> is negative.
	 *              </li>
	 *	            <li>
	 *	             If <i>m<n</i>.
	 *	            </li>
	 *	            <li>
	 *	             If the minuend is smaller than the subtrahend.
	 *	            </li>
	 *             </ul>
	 */
	static void sub
	( register uint32_t *c ,
	  register const uint32_t *a , register int m ,
	  register const uint32_t *b , register int n ) {

		register int carry;
		register int64_t tmp;
		register int j;

		// Keeps track of the carry of the 'i'th coefficient
		carry = 0;

		// Compute the first 'n' coefficients of the difference.
		// Note that 'm>=n'.
		for ( j = 0 ; j < n ; j++ , a++ , b++ , c++ ) {

			tmp = (int64_t)(*a)-(*b)+(int64_t)carry;

			if ( tmp < 0 ) {
				carry = -1;
				tmp += (((int64_t)1)<<32);
			} else {
				carry = 0;
			}

			*c = (uint32_t)(tmp&0xFFFFFFFF);
		}

		// Compute the last 'm-n' coefficients by going with and
		// updating the carry frequently.
		for ( j = n ; j < m ; j++ , a++ , c++ ) {

			tmp = (int64_t)(*a)+(int64_t)carry;

			if ( tmp < 0 ) {
				carry = -1;
				tmp += (((int64_t)1)<<32);
			} else {
				carry = 0;
			}

			*c = (uint32_t)(tmp&0xFFFFFFFF);
		}

		// If the 'carry' is non-zero then the subtrahend is smaller than the
		// minuend in which case the result is not correct.
	}

	/**
	 * @brief
	 *            Low-level function that forms the product of a
	 *            multi-precision integer by an unsigned 32-bit integer.
	 *
	 * @details
	 *            The integer
	 *            \f[
	 *             \sum_{i=0}^{n-1}a[i]\cdot 2^{i\cdot 32}
	 *            \f]
	 *            is multiplied by <i>b</i> and stored in <i>c</i> such that
	 *            \f[
	 *             carry\cdot 2^{n\cdot 32}+
	 *             \sum_{i=0}^{n-1}c[i]\cdot 2^{i\cdot 32}
	 *            \f]
	 *            is the product. The value <i>carry</i> will be returned by
	 *            the function.
	 *            <br><br>
	 *            <b>Note:</b> The arrays <i>c</i> and <i>a</i> are allowed
	 *            to be equal.
	 *
	 * @param c
	 *            Will contain the first <i>n</i> coefficients of the product
	 *            after the function has been called.
	 *
	 * @param a
	 *            Contains <i>n</i> valid unsigned 32-bit integers which
	 *            specifies the first factor.
	 *
	 * @param b
	 *            A single unsigned 32-bit integer which specifies the second
	 *            factor.
	 *
	 * @param n
	 *            Specifies the number of valid unsigned 32-bit integers that
	 *            are contained in <i>a</i>. Furthermore, <i>c</i> must be
	 *            allocated such that it can store <i>n</i> coefficients.
	 *
	 * @returns
	 *            Returns the <i>n</i>th coefficient of the product which
	 *            will not be stored in <i>c</i>.
	 *
	 * @warning
	 *            If <i>a</i> does not contain <i>n</i> valid
	 *            <code>uint32_t</code>s or if <i>c</i> cannot was not
	 *            allocated such that it can store <i>n</i> unsigned 32-bit
	 *            integers the behavior of the function is undocumented.
	 *            Furthermore, the behavior is undocumented if <i>a</i> and
	 *            <i>c</i> overlap in their respective first 'n' entries in
	 *            case they are not equal.
	 */
	static uint32_t mul
	( register uint32_t *c , register uint32_t *a , register uint32_t b ,
	  register int n ) {

		register uint32_t carry = 0;
		register uint64_t e;
		register int i;

		// Compute the first 'n' coefficients of the product by going
		// with and updating the carry
		for ( i = 0 ; i < n ; i++ , a++ , c++ ) {
			e = (uint64_t)(*a) * (uint64_t)b+carry;
			*c = (uint32_t)(e&0xFFFFFFFF);
			carry = (e>>32);
		}

		// Return the 'n'th coefficient.
		return carry;
	}

	/**
	 * @brief
	 *            Low-level function that forms the radix-\f$2^{32}\f$ product
	 *            of two positive big integers.
	 *
	 * @details
	 *            The method computes the product of the two factors
	 *            \f[
	 *             \sum_{i=0}^{m-1}a[i]\cdot 2^{i\cdot 32}
	 *            \f]
	 *            and
	 *            \f[
	 *             \sum_{i=0}^{n-1}b[i]\cdot 2^{i\cdot 32}
	 *            \f]
	 *            and stores the result in <i>c</i> such that
	 *            \f[
	 *             \sum_{i=0}^{m+n}c[i]\cdot 2^{i\cdot 32}
	 *            \f]
	 *            is the product.
	 *            <br><br>
	 *            <b>Note:</b> The array <i>c</i>, <i>a</i>, or <i>b</i>
	 *            are <b>NOT ALLOWED</b> to be pairwise equal or to
	 *            overlap/cross their first <i>n+m</i>, <i>m</i>,
	 *            or <i>n</i> elements, respectively.
	 *
	 * @param c
	 *            After calling the method, the first <i>n+m+1</i> entries of
	 *            the array contains the product as described in the details.
	 *
	 * @param a
	 *            Contains <i>m</i> valid <code>uint32_t</code>s which
	 *            specifies the first factor.
	 *
	 * @param m
	 *            The number of <code>uint32_t</code>s in <i>a</i> that
	 *            specify the first factor.
	 *
	 * @param a
	 *            Contains <i>m</i> valid <code>uint32_t</code>s which
	 *            specifies the second factor.
	 *
	 * @param n
	 *            The number of <code>uint32_t</code>s in <i>b</i> that
	 *            specify the second factor.
	 *
	 * @warning
	 *            If not sufficient memory could be provided the method
	 *            prints an error message to <code>stderr</code> and
	 *            exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            In the following cases the behavior of the function is
	 *            undocumented.
	 *            <ul>
	 *             <li>
	 *             If <i>a</i> or <i>b</i> do not contain at least <i>m</i>
	 *             or <i>n</i> valid <code>uint32_t</code>s.
	 *             </li>
	 *             <li>
	 *              If <i>c</i> was not allocated successfully to store at
	 *              least <i>m+n+1</i> unsigned 32-bit integers.
	 *             </li>
	 *             <li>
	 *              If <i>c</i>, <i>a</i>, or <i>b</i> overlap/cross their
	 *              first <i>n+m</i>, <i>m</i>, or <i>n</i> elements,
	 *              respectively, in particula, if they are not pairwise
	 *              distinct.
	 *             </li>
	 *            </ul>
	 */
	static void mul
	( register uint32_t *c ,
	  register uint32_t *a , register int m ,
	  register uint32_t *b , register int n ) {

		register uint32_t *tmp;
		register int i;

		// Allocate a temporary array that can hold 'm+1' unsigned 32-bit
		// integers
		if ( (tmp=(uint32_t*)malloc((m+1)*sizeof(uint32_t))) == NULL ) {
			cerr << "BigInteger::mul: Out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// Initialize the output array 'c' by' zeros
		memset(c,0,(m+n+1)*sizeof(uint32_t));

		// The following is an implementation of the traditional multiplication
		// scheme we learned in school with intensively use of pointer
		// arithmetic.
		b += n-1;
		c += n-1;
		for ( i = n-1 ; i >= 0 ; i-- , b-- , c-- ) {
			tmp[m] = mul(tmp,a,*b,m);
			c[m+1] += add(c,c,m+1,tmp,m+1);
		}

		// Free temporary memory
		free(tmp);
	}

	/*
	 * @param Q array allocated <code>m+1</code> integers of 32 bits. This array will
	 *          contain the radix-\f$2^{32}\f$ quotient of the input <code>U</code> by
	 *          the input <code>V</code> when terminated. More precisely
	 *          \f[
	 *              \sum_{i=0}^m Q[i]\cdot 2^{i\cdot 32}=
	 *              \lfloor
	 *                  \frac{ \sum_{i=0}^{m+n} U[i] \cdot 2^{i\cdot 32} }
	 *                       { \sum_{i=0}^{n-1} V[i] \cdot 2^{i\cdot 32} }
	 *              \rfloor.
	 *          \f]
	 *
	 *
	 * @param U on input: the numerator the radix-\f$2^{32}\f$ quotient will be formed;
	 *          on output: integral remainder times the returned value;
	 *          this array contains <code>n+m+1</code> numbers of length 32 bit and represent
	 *          the following number on input as well as on output, respectively:
	 *          \f[
	 *             \sum_{i=0}^{n+m} U[i] \cdot 2^{i\cdot 32}
	 *          \f].
	 *          On input, it is required that the most significant 32-bit word is zero, i.e.
	 *          <code>U[n+m]==0</code>.
	 *
	 * @param V the denominator the radix-\f$2^{32}\f$ quotient is formed;
	 *          this array contains <code>n</code> number of 32-bit length and represents
	 *          the following number:
	 *          \f[
	 *             \sum_{i=0}^{n-1} V[i] \cdot 2^{i\cdot 32}.
	 *          \f]
	 *          Furthermore, it is required that the most significant 32-bit is non-zero, i.e.
	 *          <code>V[n-1]!=0</code>.
	 *
	 * @param m <code> >= 0</code> the quotient <code>Q</code> will consist of at most
	 *          <code>m+1</code> words of 32-bit length
	 *
	 * @param n <code> >=2</code> the size of 32-bit words the denominator consists of;
	 *
	 * @return
	 *        32-bit integer <code>d</code> such that <code>U/d</code> is the integral
	 *          remainder of the input <code>U</code> by <code>V</code>.
	 */

	/**
	 * @brief
	 *            Low-level function that forms the radix-\f$2^{32}\f$
	 *            quotient of two positive integers.
	 *
	 * @details
	 *            The method computes the first <i>m</i> elements of <i>q</i>
	 *            such that
	 *            \f[
	 *             \sum_{i=0}^mq[i]\cdot 2^{i\cdot 32}=
	 *             \left\lfloor
	 *              \left(
	 *               \sum_{i=0}^{m+n}a[i]\cdot 2^{i\cdot 32}
	 *              \right) /
	 *              \left(
	 *               \sum_{i=0}^{n-1}b[i]\cdot 2^{i\cdot 32}
	 *              \right)
	 *             \right\rfloor
	 *            \f]
	 *            We assume from the input that \f$b[n-1]\neq0\f$ and that
	 *            \f$n\geq 2\f$ which can be asserted through binary left
	 *            shifts by 32 bits. Moreover, we assume that \f$a[m+n]=0\f$
	 *            which, again, is no real restriction for we can simply
	 *            append a zero at the end of a copy of <i>a</i>.
	 *            <br><br>
	 *            <b>Note:</b> The arrays <i>c</i>, <i>a</i>, or <i>b</i>
	 *            are <b>NOT ALLOWED</b> to be pairwise equal or to
	 *            overlap/cross their first <i>m+1</i>, <i>n+m+1</i>,
	 *            or <i>n</i> elements, respectively. Furthermore, note that
	 *            the arrays <i>a</i> and <i>b</i> will change during the
	 *            computation while the output will be the quotient
	 *            corresponding to the array <i>a</i> and <i>b</i> on
	 *            input.
 	 *            <br><br>
	 *            The algorithm was adopted from the following well-known
	 *            source in the literature:
	 *            <ul>
	 *            <li>
	 *             <b>Knuth D. E. (1981)</b>. <i>The Art of Computer
	 *             Programming, Volume II: Seminumerical Algorithms, 2nd
	 *             Edition</i>. Addison-Wesley.
	 *            </li>
	 *            </ul>
	 *
	 * @param q
	 *            Will store the radix-\f$2^{32}\f$ quotient after this
	 *            function has been called.
	 *
	 * @param a
	 *            Contains the first <i>m+n+1</i> unsigned 32-bit integers
	 *            of the numerator whose leading coefficient must be zero,
	 *            i.e. \f$b[m+n]=0\f$.
	 *
	 * @param b
	 *            Contains the first <i>n</i> unsigned 32-bit integers of the
	 *            denominator whose leading coefficient must be non-zero, i.e.
	 *            \f$b[n-1]\neq 0\f$.
	 *
	 * return
	 *           32-bit integer <i>d</i> such that <i>a/d</i> is the integral
	 *           remainder of the input <i>a</i> divided by <i>b</i>.
	 *
	 * @warning
	 *           If the assumptions as described in the details are not
	 *           satisfied, the behavior of the function is undocumented.
	 *
	 * @warning
	 *           If not sufficient memory could be provided the function
	 *           prints an error message to <code>stderr</code> and exits
	 *           with status 'EXIT_FAILURE'.
	 */
	static uint32_t div
	( uint32_t *q , uint32_t *a , uint32_t *b , int m , int n ) {

		// Note, unlike to the notation in Knuth (1981) our integer's most
		// significant coefficients are at the end of the corresponding arrays
		// and not at the beginning.

		uint32_t d , dq;
		uint32_t *db;
		uint64_t tmp , cmpl , cmpr;
		int k , j , i;

		if ( (db=(uint32_t*)malloc(sizeof(uint32_t)*(n+1))) == NULL ) {
			cerr << __FILE__ << "(" << __LINE__ << "): Out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// D1. [Normalize]
		d = (uint32_t)((((uint64_t)0xFFFFFFFF)+1) / (((uint64_t)b[n-1])+1));
		mul(a,a,d,n+m+1);
		mul(b,b,d,n);

		for ( k = 0 ; k <= m ; k++ ) {

			j = n+m-k;

			cmpr = a[j];
			cmpr <<= 32;
			cmpr += (uint64_t)a[j-1];


			// ***************************************************************
			// *************** BEGIN: D3. [Calculate 'dq'.] ******************
			// ***************************************************************

			if ( a[j] == b[n-1] ) {
				dq = 0xFFFFFFFF;
			} else {
				dq = (uint32_t)(cmpr / (uint64_t)b[n-1]);
			}

			// The loop eliminate most cases in where 'dq' is one too large
			// and all cases in where 'dq' is two too large.
			for( i = 0 ; ; ) {

				cmpl = b[n-2];
				cmpl *= (uint64_t)dq;

				tmp = b[n-1];
				tmp *= (uint64_t)dq;

				cmpr -= tmp;
				if ( cmpr > 0xFFFFFFFF ) {
					break;
				}

				cmpr <<= 32;
				cmpr += a[j-2];

				if ( cmpl <= cmpr ) {
					break;
				}

				--dq;

				if ( ++i == 2 ) {
					break;
				}

				cmpr = a[j];
				cmpr <<= 32;
				cmpr += (uint64_t)a[j-1];
			}

			// ***************************************************************
			// *************** END: D3. [Calculate 'dq'.] ********************
			// ***************************************************************

			// D4. [Multiply and subtract.]

			// Unlike Knuth (1981) we determine the correct 'dq' in a somewhat
			// more straightforward way by ensuring that 'dq*b' is smaller than
			// or equals than the integer given by the upper 'm-k' coefficients
			// of 'a'.
			db[n] = mul(db,b,dq,n);
			if ( compare(a+m-k,n+1,db,n+1) < 0 ) {
				--dq;
				db[n] = mul(db,b,dq,n);
			}
			sub(a+m-k,a+m-k,n+1,db,n+1);

			q[m-k] = dq;
		}

		free(db);

		return d;
	}

	/**
	 * @brief
	 *            Normalizes the member <code>\link size\endlink</code> to
	 *            be of minimal absolute value such that the same integer
	 *            is represented by this instances.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	void BigInteger::normalize() {

		// Index of most relevant coefficient.
		int l = getNumWords()-1;

		// Determine index of most relevant non-zero coefficient.
		for (  ; l > 0  && this->data[l] == 0 ; l-- );

		// The length is one larger
		++l;

		// If the sign is negative this is encoded by a negative
		// number of relevant words.
		if ( this->sign() < 0 ) {
			l = -l;
		}

		// Update the number of relevant coefficients.
		this->size = l;
	}

	/**
	 * @brief
	 *            Creates an instance of an big integer intialized by
	 *            the specified value.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	BigInteger::BigInteger( int64_t v ) {

		// Check whether specified integer is in the range
		// -4294967295...4294967295.
		if ( (uint64_t)(v<0?-v:v) > (uint64_t)0xFFFFFFFF ) {
			cerr << "BigInteger: "
				 << "int64_t too large." << endl;
			exit(EXIT_FAILURE);
		}

		// Allocate memory for one uint32_t
		this->data = (uint32_t*)malloc( sizeof(uint32_t) );
		if ( this->data == NULL ) {
			cerr << "BigInteger: Out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// Let the only coefficient be the absolute value
		// of the specified value.
		this->data[0] = (uint32_t)(v<0?-v:v);
		// The length is -1 if this integer is negative and
		// 1 if it is positive.
		this->size = v<0?-1:1;
		// 'data' can hold one value.
		this->capacity = 1;

	}

	/**
	 * @brief
	 *            Copy constructor.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	BigInteger::BigInteger( const BigInteger & z ) {

		// Initialize this integer by zero and ...
		this->data = (uint32_t*)malloc( sizeof(uint32_t) );
		if ( this->data == NULL ) {
			cerr << "BigInteger: "
				 << "Out of memory." << endl;
			exit(EXIT_FAILURE);
		}
		this->data[0] = 0;
		this->size = 1;
		this->capacity = 1;

		// ... use the assgnment operator for to copy
		// the integer.
		*this = z;
	}

	/**
	 * @brief
	 *            Destructor which frees all the memory that is used
	 *            to represent the integer.
	 */
	BigInteger::~BigInteger() {
		free(this->data);
	}

	/**
	 * @brief
	 *            Returns the number of bits needed to represent the absolute
	 *            value of this integer.
	 *
	 * @details
	 *            More precisely, if the absolute value of this integer is
	 *            <i>a</i> then this function computes the minimal <i>l</i>
	 *            such that \f$2^l>a\f$.
	 *
	 * @return
	 *            Number of bits to represent the absolute value of this
	 *            integer.
	 */
	int BigInteger::numBits() const {

		register int l;
		register uint32_t lw;

		l = (this->getNumWords()-1) * 32;
		lw = this->data[getNumWords()-1];

		while ( lw != 0 ) {
			++l;
			lw >>= 1;
		}

		return l;
	}

	/**
	 * @brief
	 *            Approximates the binary logarithm of this big integer.
	 *
	 * @return
	 *            Approximation of the binary logarithm of this big
	 *            integer.
	 */
	float BigInteger::log2() const {

		// Integral part of the log
		int n = numBits()-1;

		BigInteger ynum , yden;
		ynum = *this;
		leftShift(yden,1,n);

		long double y = expand(ynum,yden);

		return (float)(n + std::log(y) / std::log(2.0L));
	}

	/**
	 * @brief
	 *            Ensures that the integer is capable to hold at least
	 *            <code>numWords</code> unsigned 32-bit integers.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	void BigInteger::ensureCapacity( int numWords ) {

		if ( numWords < 0 ) {
			cerr << "BigInteger::ensureCapactiy: Capacity must be positive."
				 << endl;
			exit(EXIT_FAILURE);
		}

		if ( numWords > this->capacity ) {

			// Reallocate and check whether reallocation was successful.
			this->data = (uint32_t*)
					realloc(this->data,numWords*sizeof(uint32_t));
			if ( this->data == NULL ) {
				cerr << "BigInteger: "
					 << "Out of memory." << endl;
				exit(EXIT_FAILURE);
			}

			// Update capacity
			this->capacity = numWords;
		}

		// Ensure that non-relevant coefficients of the integer
		// are overwritten by zeros.
		memset(this->data+getNumWords(),0,
				(this->capacity-getNumWords())*sizeof(uint32_t));
	}

	/**
	 * @brief
	 *            Assign the big integer by the integer represented
	 *            by <i>a</i>.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	BigInteger &BigInteger::operator=( const BigInteger & integer ) {

		// On self-assignment we have nothing to do
		if ( this == &integer ) {
			return *this;
		}

		// Ensure this big integer can hold the relevant coefficients
		// of the integer that is copied.
		ensureCapacity(integer.getNumWords());

		// Adopt the size.
		this->size = integer.size;

		// Adopt the coefficients.
		memcpy
			(this->data,integer.data,integer.getNumWords()*sizeof(uint32_t));

		return *this;
	}

	/**
	 * @brief
	 *            Swaps the (private) member variables of this instance
	 *            with the member variables of <i>a</i>.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	void BigInteger::swap( BigInteger & a ) {

		uint32_t *data;
		int size , capacity;

		data = this->data;
		size = this->size;
		capacity = this->capacity;

		this->data = a.data;
		this->size = a.size;
		this->capacity = a.capacity;

		a.data = data;
		a.size = size;
		a.capacity = capacity;
	}

	/**
	 * @brief
	 *            Converts this integer to a long double value.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	long double BigInteger::toFloat() const {

		long double v = 0.0L;

		// First, let 'v' be the absolute value of this integer by letting
		// v = data[0]+data[1]*4294967296+...+data[i]*4294967296^i+...
		int l = getNumWords();
		for ( int i = l-1 ; i >= 0 ; i-- ) {
			v *= 4294967296.0L;
			v += (long double)this->data[i];
		}

		// Adopt the sign.
		if ( sign() < 0 ) {
			v = -v;
		}

		return v;
	}

	/**
	 * @brief
	 *             Initializes this object to represent a random integer
	 *             being smaller than the absolute value of the specified
	 *             bound.
	 *
	 * @param bound
	 *             Specifies the bound such that this object is being
	 *             assigned with an integer smaller than the absolute
	 *             value of <code>bound</code>.
	 *
	 * @param tryRandom
	 *             Indicates whether the method is advised to use
	 *             a cryptographic random generator.
	 *
	 * @warning
	 *             If not sufficient memory could be allocated, then
	 *             an error message is printed to <code>stderr</code>
	 *             and an exit with status 'EXIT_FAILURE' is caused.
	 */
	void BigInteger::random
	( const BigInteger & bound , bool tryRandom ) {

		if ( bound.isZero() ) {
			cerr << "BigInteger::random: there is no positive number being "
				 << "smaller than 0." << endl;
			exit(EXIT_FAILURE);
		}

		if ( this == &bound ) {
			BigInteger tbound(bound);
			random(tbound,tryRandom);
			return;
		}

		int size = bound.getNumWords();

		this->ensureCapacity(size);

		for ( int i = 0 ; i < size ; i++ ) {
			this->data[i] = MathTools::rand32(tryRandom);
		}
		this->size = size;
		this->normalize();

		rem(*this,*this,bound);
	}

	/**
	 * @brief
	 *             Returns the number of bytes needed to represent this
	 *             integer.
	 *
	 * @details
	 *             Let this integer be equals
	 *             \f[
	 *              \pm\sum_{j=0}^bc_j2^j
	 *             \f]
	 *             where \f$c_j\in\{0,1\}\f$. Then this function returns
	 *             \f$\lceil (b+1)/8 \rceil\f$ which is the number of
	 *             bytes needed to represent this integer including a
	 *             leading bit encoding its sign.
	 *
	 * @return
	 *             The number of bytes needed to represent this integer.
	 */
	int BigInteger::getSizeInBytes() const {

		int tmp = numBits()+1;
		return tmp/8+(tmp%8?1:0);
	}

	/**
	 * @brief
	 *             Exports this big integer to byte data.
	 *
	 * @details
	 *             This method write bytes to <code>array</code> such
	 *             that the value
	 *             \f[
	 * (array[s-1]\&127)\cdot 256^{s-1}+\sum_{i=0}^{s-2} array[j]\cdot 256^j
	 *             \f]
	 *             equals the absolute value of this integer and
	 *             <i>s=</i>\link getSizeInBytes()\endlink. If this
	 *             is negative, then the leading bit of
	 *             <code>array[s-1]</code> is set to one.
	 *
	 *             To obtain a BigInteger object back from the exported
	 *             data, we may use the \link fromBytes()\endlink
	 *             method.
	 *
	 * @param array
	 *             The byte array into which the data of this integer
	 *             is exported.
	 *
	 * @see fromBytes()
	 * @see getSizeInBytes()
	 *
	 * @warning
	 *             If <code>array</code> is not capable of holding at
	 *             least \link getSizeInBytes()\endlink bytes, then
	 *             calling this method results into undocumented
	 *             behavior.
	 */
	void BigInteger::toBytes( uint8_t *array ) const {

		int sizeInBytes = getSizeInBytes();

		array[sizeInBytes-1] = 0;

		int k = 0;
		for ( int i = 0 ; i < getNumWords() ; i++ ) {

			uint32_t v = this->data[i];

			for ( int j = 0 ; j < 4 && k < sizeInBytes ; j++ , k++ ) {
				array[k] = v & 0xFF;
				v >>= 8;
			}
		}

		if ( this->size < 0 ) {
			array[sizeInBytes-1] += 128;
		}
	}

	/**
	 * @brief
	 *             Initializes this big integer by the data contained
	 *             in the specified array.
	 *
	 * @details
	 *             This object will be initialized such that its
	 *             absolute value will be equals
	 *             \f[
	 * (array[s-1]\&127)\cdot 256^{s-1}+\sum_{i=0}^{s-2} array[j]\cdot 256^j
	 *             \f]
	 *             where <i>s=</i><code>sizeInBytes</code>. If the leading
	 *             bit of <code>array[s-1]</code> is set, then the integer
	 *             will be negative.
	 *
	 *             This method can be considered as the inverse of
	 *             the \link toBytes()\endlink method.
	 *
	 * @param array
	 *             The array containing the data from which this big
	 *             integer is initialized.
	 *
	 * @param sizeInBytes
	 *             The number of valid bytes stored in
	 *             <code>array</code>.
	 *
	 * @see toBytes()
	 *
	 * @warning
	 *             If <code>array</code> does not contain
	 *             <code>sizeInBytes</code> valid bytes, then calling
	 *             this method runs into undocumented behavior.
	 *
	 * @warning
	 *             If <code>sizeInBytes</code> is negative, an error
	 *             message will be printed to <code>stderr</code>
	 *             and the program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message will be printed to <code>stderr</code>
	 *             and the program exits with status 'EXIT_FAILURE'.
	 */
	void BigInteger::fromBytes( const uint8_t *array , int sizeInBytes ) {

		if ( sizeInBytes < 0 ) {
			cerr << "BigInteger::fromBytes: number of bytes "
				 << "must be positive." << endl;
			exit(EXIT_FAILURE);
		}

		clear();

		if ( sizeInBytes != 0 ) {

			int numWords = sizeInBytes/4 + (sizeInBytes%4?1:0);

			ensureCapacity(numWords);
			this->size = numWords;

			for ( int i = 0 ; i < numWords ; i++ ) {
				this->data[i] = 0;
				for ( int j = 3 ; j >= 0 ; j-- ) {
					if ( i*4+j < sizeInBytes ) {
						this->data[i] <<= 8;
						if ( i*4+j == sizeInBytes-1 ) {
							this->data[i] += array[i*4+j]&127;
						} else {
							this->data[i] += array[i*4+j];
						}
					}
				}
			}

			normalize();

			if ( array[sizeInBytes-1] & 128 ) {
				negate(*this,*this);
			}
		}
	}

	/**
	 * @brief
	 *             Initializes this big integer by a string in
	 *             decimal representation.
	 *
	 * @details
	 *             If the specified string <code>str</code> contains
	 *             characters that are no digits (except a possible
	 *             leading '-' or '+'), then the integer will be
	 *             initialized by the first digits only and the
	 *             function returns <code>false</code>; otherwise
	 *             if <code>str</code> represents a valid decimal
	 *             integer, the function returns <code>true</code>.
	 *
	 * @param str
	 *             The string with which this big integer is initialized.
	 *
	 * @return
	 *             <code>true</code> if <code>str</code> represents
	 *             a valid decimal integer; otherwise, the function
	 *             returns <code>false</code>.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message will be printed to <code>stderr</code>
	 *             and the program exits with status 'EXIT_FAILURE'.
	 */
	bool BigInteger::fromString( const string & str ) {

		clear();

		int n = (int)str.length();
		const char *c_str = str.c_str();
		int s = 1;

		for ( int i = 0 ; i < n ; i++ ) {

			if ( i == 0 && (c_str[0] == '-' || c_str[0] == '+') ) {

				if ( c_str[0] == '-' ) {
					s = -1;
				}

				if ( n == 1 ) {
					return false;
				}

			} else {

				if ( !isdigit(c_str[i]) ) {
					return false;
				}

				mul(*this,*this,BigInteger(10));
				add(*this,*this,BigInteger( (int64_t)(c_str[i]-'0') ) );
			}

		}

		if ( s < 0 ) {
			negate(*this,*this);
		}

		return true;
	}

	/**
	 * @brief
	 *            Compares the absolute values of two integers.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	int BigInteger::absCompare
	( const BigInteger & a , const BigInteger & b ) {

		int m , n;
		m = a.getNumWords();
		n = b.getNumWords();

		// If the number of relevant coefficients already allow to decide
		// which integer is larger and smaller, return correspondingly.
		if ( m > n ) {
			return 1;
		} else if ( m < n ) {
			return -1;
		}

		// Here, 'a' and 'b' both have the same number of relevant
		// coefficients. Thus, return result of corresponding low-level
		// function.
		return thimble::compare(a.data,m,b.data,n);
	}

	/**
	 * @brief
	 *            Compares two signed big integer.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	int BigInteger::compare
	( const BigInteger & a , const BigInteger & b ) {

		// The first steps are to check whether the signs of both integers
		// are different which allows to return immediately.

		if ( a.sign() > 0 && b.sign() < 0 ) {
			return 1;
		}

		if ( a.sign() < 0 && b.sign() > 0 ) {
			return -1;
		}

		if ( a.sign() == 0 && b.sign() == 0 ) {
			return 0;
		}

		if ( a.sign() == 0 ) {
			return -b.sign();
		}

		if ( b.sign() == 0 ) {
			return a.sign();
		}


		// Now, the signs are equal and ..
		int s = a.sign();

		// ... the result of the absolute comparison is signed
		// correspondingly.
		return s * absCompare(a,b);
	}

	/**
	 * @brief
	 *            Performs a binary left shift of a big integer by
	 *            <i>n</i> bits.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	void BigInteger::leftShift
	( BigInteger & a , const BigInteger & b , int n ) {

		// If output and input are of same reference, it is necessary to
		// compute the left shift of a copy.
		if ( &a == &b ) {
			BigInteger tb = b;
			leftShift(a,tb,n);
			return;
		}

		// If the specified bit positions is negative, the left shift
		// is interpreted as a right shift.
		if ( n < 0 ) {
			rightShift(a,b,-n);
			return;
		}

		// Backup of the sign.
		int s = b.sign();

		// Left shift by 'n=k*32+m' bits have to be performed.
		int k , m;
		k = n/32L;
		m = n%32L;

		// Temporarily, let 'a' be of as many zero words which suffice to
		// hold the result; we will normalize later.
		a.clear();
		a.ensureCapacity(k+b.getNumWords()+1);
		a.size =  k+b.getNumWords()+1;

		// The first 'k' coefficients, which are 32-bit words, of 'a' will be
		// zero because 'k' is the proportion of whole 32 shift bit positions;
		// thus, it suffices to compute the coefficients of the output which
		// are at index >='k'.
		for ( int i = b.getNumWords() - 1 ; i >= 0 ; i-- ) {

			// The first '32-m' bits of the current coefficient will be
			// stored in 'u'...
			uint64_t u = (uint64_t)b.data[i];
			u <<= m;

			// ... and be combined with the last 'm' bits.
			a.data[k+i+1] |= (u>>32);
			a.data[k+i] |= (u&0xFFFFFFFF);
		}

		// There is no harm in normalizing the shift.
		a.normalize();

		// Adopt the backuped sign.
		if ( s < 0 ) {
			a.size = -a.size;
		}
	}

	/**
	 * @brief
	 *            Performs a binary right shift of a big integer by
	 *            <i>n</i> bits.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	void BigInteger::rightShift
	( BigInteger & a , const BigInteger & b , int n ) {

		// If output and input are of same reference, it is necessary to
		// compute the right shift of a copy.
		if ( &a == &b ) {
			BigInteger tb = b;
			rightShift(a,tb,n);
			return;
		}

		// If the specified bit positions is negative, the right shift
		// is interpreted as a left shift.
		if ( n < 0 ) {
			leftShift(a,b,-n);
			return;
		}

		// Backup of the sign.
		int s = b.sign();

		// Length of unsigned 32-bit integers needed to represent the absolute
		// value of 'b'
		int blen = b.getNumWords();

		// Right shift by 'n=k*32+m' bits have to be performed.
		int k , m;
		k = n/32L;
		m = n%32L;

		// Temporarily, let 'a' be of as many zero words which suffice to
		// hold the result; we will normalize later.
		a.clear();
		a.ensureCapacity(blen-k);
		a.size = blen-k;

		// Adopt the right-shifted bits by computing the 'k' coefficients
		// of the result.
		for ( int i = k ; i < blen ; i++ ) {
			a.data[i-k] = (b.data[i]>>m);
			if ( i+1 < blen ) {
				a.data[i-k] += (b.data[i+1]&((1<<m)-1))<<(32-m);
			}
		}

		// There is no harm in normalizing the shift.
		a.normalize();

		// Adopt the backuped sign.
		if ( s < 0 ) {
			negate(a,a);
		}
	}

	/**
	 * @brief
	 *            Computes the sum of two big integers.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	void BigInteger::add
	( BigInteger & c , const BigInteger & a , const BigInteger & b ) {


		// Trivial cases.
		if ( b.sign() == 0 ) {
			c = a;
			return;
		} else if ( a.sign() == 0 ) {
			c = b;
			return;
		}

		if ( (a.sign() >= 0 && b.sign() >= 0) ||
			 (a.sign() < 0 && b.sign() < 0) ) {

			// ***************************************************************
			//  If the signs are equal, we wrap around the low-level function
			// 'add(uint32_t*,uint32_t*,int,uint32_t*,int)' to compute the
			// result.
			// ***************************************************************

			// Sign of the output
			int s = a.sign();

			int m , n;
			m = a.getNumWords();
			n = b.getNumWords();

			// Ensure, the capacity of the output is large enough to hold
			// the sum.
			int l;
			l = std::max(m,n)+1;
			c.ensureCapacity(l);

			// Temporarily, let the number of relevant words be such that
			// it can hold the sum, even on a possible carry of the last
			// coefficient --- we will normalize later.
			c.size = l;

			// Use the low-level function to compute the sum. The 'l-1'th
			// coefficient can be both, zero and non-zero. Therefore, ...
			c.data[l-1] = thimble::add(c.data,a.data,m,b.data,n);

			// ..., we normalize the sum.
			c.normalize();

			// If the signs of the summands were negative, let the sign of
			// the sum be negative.
			if ( s < 0 ) {
				c.size = -c.size;
			}

		} else {

			// ***************************************************************
			// If the signs are not equal, the addition essentially is a
			// subtraction; thus we wrap around the low-level function
			// 'sub(uint32_t*,uint32_t*,int,uint32_t*,int) to compute
			// the result.
			// ***************************************************************

			// For the low-level subtraction function, we need to known which
			// if the summands is larger in its absolute size.
			int ac = absCompare(a,b);

			// If both summands are equal in its absolute size, the result will
			// be zero, because the sum is essentially a difference for the
			// signs of the input differ. ...
			if ( ac != 0 ) {

				// ... . If the absolute values differ, we have to perform the
				// following.

				int m , n;
				m = a.getNumWords();
				n = b.getNumWords();

				// The sign of the output
				int s = a.sign()*ac;

				// Temporarily, let the first 'max(m,n)' words of 'c' be
				// relevant. We will normalize later.
				int l;
				l = max(m,n);
				c.ensureCapacity(l);
				c.size = l;

				// We use the low-level subtraction function to compute the
				// difference which requires the minuend to be of larger
				// absolute value than the subtrahend. Therefore, we must
				// consider the two possible cases. The sign, which was
				// computed above, will guarantee the correctness of
				// the output integer 'c'.
				if ( ac > 0 ) {
					thimble::sub(c.data,a.data,m,b.data,n);
				} else {
					thimble::sub(c.data,b.data,n,a.data,m);
				}

				// Leading coefficients may have become zero. Therefore we
				// must normalize.
				c.normalize();

				// Adopt the sign.
				if ( s < 0 ) {
					negate(c,c);
				}

			}
		}
	}

	/**
	 * @brief
	 *            Computes the difference between two big integers.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	void BigInteger::sub
	( BigInteger & c , const BigInteger & a , const BigInteger & b ) {

		// Trivial cases.
		if ( a.sign() == 0 ) {
			negate(c,b);
			return;
		} else if ( b.sign() == 0 ) {
			c = a;
			return;
		}

		if ( (a.sign() > 0 && b.sign() > 0) ||
			 (a.sign() < 0 && b.sign() < 0) ) {

			// ***************************************************************
			//  If the signs are equal, we wrap around the low-level function
			// 'sub(uint32_t*,uint32_t*,int,uint32_t*,int)' to compute the
			// result.
			// ***************************************************************

			int ac = absCompare(a,b);

			// If the absolute values differ, we perform the following.
			if ( ac != 0 ) {

				int m , n;
				m = a.getNumWords();
				n = b.getNumWords();

				// The sign of the "difference".
				int s = ac * a.sign();

				// Temporarily, let the first 'max(m,n)' words of 'c' be
				// relevant. We will normalize later.
				int l = max(m,n);
				c.ensureCapacity(l);
				c.size = l;

				// We use the low-level subtraction function to compute the
				// difference which requires the minuend to be of larger
				// absolute value than the subtrahend. Therefore, we must
				// consider the two possible cases. The sign, which was
				// computed above, will guarantee the correctness of
				// the output integer 'c'.
				if ( ac > 0 ) {
					thimble::sub(c.data,a.data,m,b.data,n);
				} else {
					thimble::sub(c.data,b.data,n,a.data,m);
				}

				// Leading coefficients may have become zero. Therefore we
				// must normalize.
				c.normalize();

				// Adopt the sign.
				if ( s < 0 ) {
					c.size = -c.size;
				}

			} else {
				// If both minuend and subtrahend are equal in its absolute size,
				// the result will be zero.
				c.clear();
			}

		} else {

			// ***************************************************************
			// If the signs are not equal, the subtraction essentially is an
			// addition; thus we wrap around the low-level function
			// 'add(uint32_t*,uint32_t*,int,uint32_t*,int) to compute
			// the result.
			// ***************************************************************

			// The sign of the output which is essentially a sum.
			int s = a.sign();

			int m , n;
			m = a.getNumWords();
			n = b.getNumWords();

			// Temporarily, let the first 'max(m,n)' words of 'c' be
			// relevant. We will normalize later to get rid of a
			// possible zero leading coefficient.
			int l = max(m,n)+1;
			c.ensureCapacity(l);
			c.size = l;

			// Use the low-level function to compute the sum. The 'l-1'th
			// coefficient can be both, zero and non-zero. Therefore, ...
			c.data[l-1] = thimble::add
				(c.data,a.data,a.getNumWords(),b.data,b.getNumWords());

			// ..., we normalize the sum.
			c.normalize();

			// Adopt the sign.
			if ( s < 0 ) {
				c.size = -c.size;
			}
		}
	}

	/**
	 * @brief
	 *            Computes the product of two big integers.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	void BigInteger::mul
	( BigInteger & c , const BigInteger & a , const BigInteger & b ) {

		// To use the low-level function
		// 'mul(uint32_t*,uint32_t*,int,uint32_t*,int)' the three arrays must
		// not be equal or overlap. Thus, we make copies if there are equal
		// references in which case we recursively call the method.
		if ( &c == &a ) {
			BigInteger ta = a;
			mul(c,ta,b);
			return;
		} else if ( &c == &b ) {
			BigInteger tb = b;
			mul(c,a,tb);
			return;
		}

		// The sign of the ouput.
		int s = a.sign()*b.sign();

		// Temporarily, let the size of the output be such that it can
		// hold the product in any case. We will normalize later.
		c.clear();
		c.ensureCapacity(a.getNumWords()+b.getNumWords()+1);
		c.size = a.getNumWords()+b.getNumWords()+1;

		// Call of low-level function.
		thimble::mul(c.data,a.data,a.getNumWords(),b.data,b.getNumWords());

		// The leading 'l-1'th coefficient can be zero.
		// Thus we must normalize.
		c.normalize();

		// Adopt the sign.
		if ( s < 0 ) {
			negate(c,c);
		}
	}

	/**
	 * @brief
	 *            Computes the quotient of two big integers.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	void BigInteger::div
	( BigInteger & c , const BigInteger & a , const BigInteger & b ) {

		// Checks for trivial cases.
		if ( b.isZero() ) {
			cerr << __FILE__ << "(" << __LINE__ << "): Division by zero." << endl;
			exit(EXIT_FAILURE);
		} if ( a.sign() > 0 && b.sign() > 0 && a.getNumWords() < b.getNumWords() ) {
			c.clear();
			return;
		} else if ( b.getNumWords() == 1 ) {
			BigInteger ta , tb;
			leftShift(ta,a,32);
			leftShift(tb,b,32);
			div(c,ta,tb);
			return;
		}

		int m , n;
		n = b.getNumWords();
		m = a.getNumWords()-n;

		// Sign of the quotient.
		int s = a.sign() * b.sign();

		// Using the low-level function
		// 'div(uint32_t*,uint32_t*,uint32_t*,int,int)' modifies the input
		// arrays. Thus, we must make of copy of the data of 'a' and 'b'
		// which requires to allocate memory.
		uint32_t *u , *v;
		u = (uint32_t*)malloc((m+n+1)*sizeof(uint32_t));
		v = (uint32_t*)malloc(n*sizeof(uint32_t));
		if ( u == NULL || v == NULL ) {
			cerr << "BigInteger::div: Out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// Copy the data of 'a' and 'b'.
		memcpy(u,a.data,(m+n)*sizeof(uint32_t));
		memcpy(v,b.data,n*sizeof(uint32_t));
		// The low-level function requires that the leading coefficient
		// of the numerator is zero.
		u[m+n] = 0;

		// Temporarily, let the size of the qutient be such that it can
		// hold the result.
		c.clear();
		c.ensureCapacity(m+1);
		c.size = m+1;

		// Call of low-level function.
		thimble::div(c.data,u,v,m,n);

		// Ensure that the leading coefficient is non-zero if the quotient
		// is non-zero.
		c.normalize();

		if ( s < 0 ) {

			// If the quotient is negative because we want the quotient
			// to be the floor of 'a/b' we check if the remainder is non-zero
			// ...
			bool isZero = true;
			for ( int i = 0 ; i < m+n ; i++ ) {
				if ( u[i] ) {
					isZero = false;
					break;
				}
			}

			//... and if it is not, ...
			if ( !isZero ) {
				// we increment the absolute value of the qutient.
				add(c,c,1);
			}

			// If the quotient is negative, negate.
			negate(c,c);
		}

		// Free memory.
		free(u);
		free(v);
	}

	/**
	 * @brief
	 *            Computes the Euclidean division of two integers.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	void BigInteger::divRem
	( BigInteger & q , BigInteger & r ,
	  const BigInteger & a , const BigInteger & b ) {

		// Check if quotient and remainder can be distinguished.
		if ( &q == &r ) {
			cerr << "BigInteger::divRem: "
				 << "Quotient must be of different reference than remainder."
				 << endl;
			exit(EXIT_FAILURE);
		}

		// After the quotient has been computed we need the input 'a' and 'b'
		// to compute the remainder. Thus, we need to ensure 'a' is of
		// different reference than 'q' and 'b' and that 'b' is of different
		// reference than 'r'.

		if ( &q == &a || &r == &a ) {
			BigInteger ta = a;
			divRem(q,r,ta,b);
			return;
		}

		if ( &q == &b ) {
			BigInteger tb = b;
			divRem(q,r,a,tb);
			return;
		}

		// q = floor(a/b)
		div(q,a,b);

		// r = a-q*b
		mul(r,q,b);
		sub(r,a,r);
	}

	/**
	 * @brief
	 *            Computes the greatest common divisor of two big
	 *            integers.
	 *
	 * @param g
	 *            On output, the greatest common divisor of <code>a</code>
	 *            and <code>b</code>.
	 *
	 * @param a
	 *            First integer.
	 *
	 * @param b
	 *            Second integer.
	 *
	 * @warning
	 *            If <code>a</code> and <code>b</code> represent zero
	 *            integer, an error message will be printed to
	 *            <code>stderr</code> and the program exits with
	 *            status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not sufficient memory could be allocated, the
	 *            method prints an error message to <code>stderr</code>
	 *            and exits with status 'EXIT_FAILURE'.
	 */
	void BigInteger::gcd
	( BigInteger & g , const BigInteger & a , const BigInteger & b ) {

		if ( a.isZero() && b.isZero() ) {
			cerr << "BigInteger::gcd: invalid arguments." << endl;
			exit(EXIT_FAILURE);
		}

        BigInteger g1 , g2 , tmp;

        g1 = b; g = a;

        divRem(tmp,g2,g,g1);

        while ( !g1.isZero() ) {
            divRem(tmp,g2,g,g1);
            g1.swap(g2);
            g.swap(g2);
        }
	}

	/**
	 * @brief
	 *            Computes the factorial of an integer as a big integer.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	BigInteger BigInteger::factorial( int n ) {

		// Check if input is reasonable.
		if ( n < 0 ) {
			cerr << "BigInteger::factor: Argument "
				 << "must be positive." << endl;
			exit(EXIT_FAILURE);
		}

		BigInteger a , b , c;

		a.data[0] = 1;
		for ( int i = 1 ; i <= n ; i++ ) {

			// b = a*i;
			c.data[0] = i;
			mul(b,a,c);

			// a<->b such that 'a<-a*i'
			a.swap(b);
		}

		return a;
	}

	/**
	 * @brief
	 *            Computes the binomial coefficient &quot;<i>n</i> over
	 *            <i>k</i>&quot; as a big integer.
	 *
	 * @details
	 *            'see BigInteger.h'
	 */
	BigInteger BigInteger::binomial( int n , int k ) {

		// Check if the inputs are reasonable.

		if ( k < 0 || n < 0 ) {
			cerr << "BigInteger::binomial: Bad arguments. Arguments "
				 << "must be positive." << endl;
			exit(EXIT_FAILURE);
		}

		if ( k > n ) {
			cerr << "BigInteger::binomial: Bad arguments. Lower term must be "
				 << "smaller than or equals the upper term." << endl;
			exit(EXIT_FAILURE);
		}

		if ( n-k < k ) {
			return binomial(n,n-k);
		}

		BigInteger num(1) , den(1) , tmp1 , tmp2;

		// num = (n+1-1)*(n+1-2)*...*(n+1-k)
		for ( int j = 1 ; j <= k ; j++ ) {
			tmp1.data[0] = n+1-j;
			mul(tmp2,num,tmp1);
			num.swap(tmp2);
		}

		// den = 1*2*...*k
		for ( int j = 1 ; j <= k ; j++ ) {
			tmp1.data[0] = j;
			mul(tmp2,den,tmp1);
			den.swap(tmp2);
		}

		// tmp1=\left({n\atop k}\right) = \prod_{j=1}^k (n+1-j)/j
		div(tmp1,num,den);

		return tmp1;
	}

	/**
	 * @brief
	 *            Expands a float from a fraction of big integers.
	 *
	 * @param num
	 *            Numerator.
	 *
	 * @param den
	 *            Denonimator.
	 *
	 * @param numBits
	 *            Number of bits with which the fractional part of
	 *            the ratio is approximated.
	 *
	 * @return
	 *            A float approximating <code>num/den</code>
	 *
	 * @warning
	 *            If <code>den</code> represents the zero number, an
	 *            error message will be printed to <code>stderr</code>
	 *            and the program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not sufficient memory could be allocated, the
	 *            method prints an error message to <code>stderr</code>
	 *            and exits with status 'EXIT_FAILURE'.
	 */
	long double BigInteger::expand
	( const BigInteger & num , const BigInteger & den ,
	  int numBits ) {

		if ( den.isZero() ) {
			cerr << "BigInteger::expand: Division by zero." << endl;
			exit(EXIT_FAILURE);
		}

		BigInteger q , r;
		divRem(q,r,num,den);

		long double R = 0.0L;

		for ( int m = 1 ; m < numBits ; m++ ) {
			leftShift(r,r,1);
			if ( compare(r,den) >= 0 ) {
				R += pow(0.5,m);
				r -= den;
			}
		}

		return R + q.toFloat();
	}

	/**
	 * @brief
	 *            Prints the big integer in decimal text representation to the
	 *            specified output stream.
	 *
	 * @details
	 *            see 'BigInteger.h'
	 */
	ostream & operator<<
			( ostream & out , const BigInteger & a ) {

		if ( a.isZero() ) {
			out << "0";
			return out;
		}

		// Allocate an array which is guaranteed to hold the number of digits
		// needed to represent the integer in decimal representation.
		unsigned char *digits = (unsigned char*)malloc(a.getNumWords() * 10);
		if ( digits == NULL ) {
			cerr << "operator<<(std::ostream&,const thimble::BigInteger: "
				 << "Out of memory" << endl;
			exit(EXIT_FAILURE);
		}

		// Keeps track of the number of decimal digits needed to represent
		// the integer.
		int n = 0;

		{
			// Collect digits by successfully computing 'r=b rem 10' until
			// quotient becomes zero.
			BigInteger b , q , r;
			BigInteger::abs(b,a);
			do {
				divRem(q,r,b,10);
				digits[n++] = '0'+(unsigned char)r.toInt();
				b.swap(q);
			} while ( !b.isZero() );
		}

		// Output '-' if integer is negative.
		if ( a.sign() < 0 ) {
			out << "-";
		}

		// Output the digits in converse order.
		for ( int i = n-1 ; i >= 0 ; i-- ) {
			out << digits[i];
		}

		// Free memory holding the digits.
		free(digits);

		return out;
	}
}

