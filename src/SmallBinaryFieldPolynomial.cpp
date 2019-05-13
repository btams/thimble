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
 * @file SmallBinaryFieldPolynomial.cpp
 *
 * @brief
 *            Implements the functions provided by
 *            'SmallBinaryFieldPolynomial.h' which is related with a class to
 *            represent and compute with polynomials having coefficients in a
 *            small binary galois field
 *            (\link thimble::SmallBinaryField\endlink).
 *
 * @author Benjamin Tams
 */

#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include <thimble/math/MathTools.h>
#include <thimble/math/numbertheory/SmallBinaryField.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *           Normalizes the degree of the polynomial
	 *
	 * @details
	 *           see 'SmallBinaryFieldPolynomial.h'
	 */
	void SmallBinaryFieldPolynomial::normalize() {

		for ( int i = this->degree ; i >= 0 ; i-- ) {
			if ( this->coefficients[i] == 0 )
				this->degree -= 1;
			else
				break;
		}
	}

	/**
	 * @brief
	 *           Assigns this polynomial with the specified polynomial.
	 *
	 * @details
	 *           see 'SmallBinaryFieldPolynomial.h'
	 */
	void SmallBinaryFieldPolynomial::assign
			( const SmallBinaryFieldPolynomial & f ) {

		if ( this == &f ) {
			return;
		}

		if ( this->gfPtr->getDefiningPolynomial().rep != f.gfPtr->getDefiningPolynomial().rep ) {
			cerr << "SmallBinaryFieldPolynomial::assign: "
				 << "related to different finite fields." << endl;
			exit(EXIT_FAILURE);
		}

		if ( f.isZero() ) {
			setZero();
			return;
		}

		ensureCapacity(f.degree+1);

		memcpy( this->coefficients , f.coefficients ,
				(f.degree+1)*sizeof(uint32_t));

		this->degree = f.degree;
	}

	/**
	 * @brief
	 *           Ensures that the polynomial has enough capacity
	 *           to hold at least <code>newCapacity</code> elements.
	 *
	 * @details
	 *           see 'SmallBinaryFieldPolynomial.h'
	 *
	 */
	void SmallBinaryFieldPolynomial::ensureCapacity( int newCapacity ) {

		// Check whether capacity must be increased
		if ( newCapacity > this->capacity ) {

			uint32_t *newCoefficients;

			// If no space was yet allocated...
			if ( this->capacity == 0 ) {

				// ... allocate space correspondingly ...
				newCoefficients = (uint32_t*)malloc
						(newCapacity*sizeof(uint32_t));
			} else {

				// ... otherwise reallocate the space.
				newCoefficients = (uint32_t*)realloc
						(this->coefficients,newCapacity*sizeof(uint32_t));
			}

			// Check whether space could be successfully allocated ...
			if ( newCoefficients == NULL ) {

				// ... if not, verbose and exit.
				cerr << "SmallBinaryFieldPolynomial::ensureCapacity: "
				     << "Out of memory" << endl;
				exit(EXIT_FAILURE);
			}

			// Initialize the newly added space with 0 bits
			memset
			( newCoefficients+this->capacity , 0 ,
			  (newCapacity-this->capacity)*sizeof(uint32_t));

			// Update the polynomial
			this->coefficients = newCoefficients;
			this->capacity = newCapacity;
		}
	}

	/**
	 * @brief
	 *           Determine whether this polynomials is equals to the
	 *           given polynomial.
	 *
	 * @details
	 *           see 'SmallBinaryFieldPolynomial.h'
	 */
	bool SmallBinaryFieldPolynomial::equals
	( const SmallBinaryFieldPolynomial & f ) const {

		// If 'f' is of the same reference as this polynomial
		// the polynomials must be equal.
		if ( this == &f ) {
			return true;
		}

		// Check whether the two polynomials are defined over the same
		// finite field.
		if ( this->gfPtr != f.gfPtr ) {
			cerr << "SmallBinaryFieldPolynomial: When two polynomials are "
				 << "compared they must have coefficients in the same field."
				 << endl;
			exit(EXIT_FAILURE);
		}

		// If the degrees differ the polynomial can not be equal.
		int n = this->degree;
		if ( n != f.degree ) {
			return false;
		}

		// Compare the relevant bytes of the coefficients ...
		int c = memcmp
				(this->coefficients,f.coefficients,(n+1)*sizeof(uint32_t));

		// ... which are equal if the result is 0.
		return c == 0;
	}

	/**
	 * @brief
	 *            Sets/Replaces the <i>i</i>-th coefficient of the
	 *            polynomial.
	 *
	 * @details
	 *            see 'SmallBinaryFieldPolynomial.h'
	 */
	void SmallBinaryFieldPolynomial::setCoeff( int i , uint32_t c ) {

		// Check whether a valid index is accessed ...
		if ( i < 0 ) {
			// ... if not, verbose and exit.
			cerr << "SmallBinaryFieldPolynomial::setCoeff: Accessed "
				 << "coeffcient must have positive index." << endl;
			exit(EXIT_FAILURE);
		}

		// Save the degree of the polynomial.
		int d = this->degree;

		if ( d < i ) {

			// If a coefficient beyond relevant indices
			// is accessed, ensure sufficient capacity, ...
			ensureCapacity(i+1);

			// ..., ensure the offset to be zero, and ...
			for ( int j = d+1 ; j < i ; j++ ) {
				this->coefficients[j] = 0;
			}

			// ..., and update the degree
			this->degree = i;
		}

		// Update the accessed coeffcient
		this->coefficients[i] = c;

		// The degree might be wrong if accessed index is
		// higher than the old degree of the polynomial
		if ( d <= i ) {
			normalize();
		}
	}

	/**
	 * @brief
	 *            Replaces this polynomial by a random polynomial
	 *            of degree smaller than <code>size</code>
	 *
	 * @details
	 *            see 'SmallBinaryFieldPolynomial.h'
	 */
	void SmallBinaryFieldPolynomial::random( int size , bool tryRandom ) {

		// Ensure the polynomial has sufficient capacity
		ensureCapacity(size);

		// If the size of the polynomial is smaller than or equal 0
		// the polynomial must be zero
		if ( size <= 0 ) {
			this->degree = -1;
			if ( this->coefficients != NULL ) {
				this->coefficients[0] = 0;
			}
			return;
		}

		// Temporarily, set the maximal possible degree
		this->degree = size-1;

		// The first 'size' coefficients are successively generated
		for ( int i = 0 ; i < size ; i++ ) {
			this->coefficients[i] = this->gfPtr->random(tryRandom);
		}

		// Normalize the degree of the polynomial
		normalize();
	}

	/**
	 * @brief
	 *            Replaces the polynomial by the polynomial that
	 *            consists of linear factors corresponding to
	 *            having the specified elements as root.
	 *
	 * @details
	 *            see 'SmallBinaryFieldPolynomial.h'
	 */
	void SmallBinaryFieldPolynomial::buildFromRoots
	( const uint32_t *a , int n ) {

		SmallBinaryFieldPolynomial tmp(getField());
		tmp.ensureCapacity(n+1);

		ensureCapacity(n+1);
		setOne();

		SmallBinaryFieldPolynomial x(this->gfPtr[0]);
		x.setX();

		for ( int i = 0 ; i < n ; i++ ) {
			x.coefficients[0] = this->gfPtr->neg(a[i]);
			mul(tmp,*this,x);
			this->swap(tmp);
		}
	}

	/**
	 * @brief
	 *            Finds all roots of the polynomial in the coefficient
	 *            field.
	 *
	 * @details
	 *            see 'SmallBinaryFieldPolynomial.h'
	 */
	int SmallBinaryFieldPolynomial::findRoots( uint32_t *x ) const {

		int d = deg();
		if ( d < 0 ) {
			cerr << "SmallBinaryFieldPolynomial::findRoots: Polynomial must "
				 << "be non-zero." << endl;
			exit(EXIT_FAILURE);
		}

		if ( d == 0 ) {
			return 0;
		}

		SmallBinaryFieldPolynomial
			f(*this) , q(*(this->gfPtr)) , r(*(this->gfPtr)) ,
			tmp(*(this->gfPtr));

		uint64_t size = this->gfPtr->getCardinality();

		int n = 0;
		uint32_t a , b;
		for ( a = 0 ; a < size ; a++ ) {
			b = f.eval(a);
			if ( b == 0 ) {

				// Remove all multiplicities of the corresponding linear
				// factor from 'f'.
				tmp.setX();
				tmp.setCoeff(0,a);
				divRem(q,r,f,tmp);
				do {
					f.swap(q);
					divRem(q,r,f,tmp);
				} while ( r.isZero() );

				x[n++] = a;

				// Check if more roots are possible.
				if ( n >= d ) {
					break;
				}
			}
		}

		return n;
	}

	/**
	 * @brief
	 *            Evaluates the polynomial at a field element.
	 *
	 * @details
	 *            see 'SmallBinaryFieldPolynomial.h'
	 */
	uint32_t SmallBinaryFieldPolynomial::eval( uint32_t x ) const {

		int d = deg();

		if ( d >= 0 ) {

			// Implementation of Horner's method if the degree
			// is larger than or equals 0
			register uint32_t y = getCoeff(d);
			register uint32_t *a = this->coefficients+d-1;
			register int i;
			for ( i = d-1 ; i >= 0 ; i-- , a-- ) {
				y = this->gfPtr->mul(y,x);
				y ^= *a;
			}

			return y;

		}

		return 0;
	}

	/**
	 * @brief
	 *            Computes the polynomial of minimal degree that
	 *            interpolates the tuples <i>(a[i],b[i])</i>
	 *            for <i>i=0,...,n-1</i>.
	 *
	 * @details
	 *            see 'SmallBinaryFieldPolynomial.h'
	 */
	void SmallBinaryFieldPolynomial::interpolate
	( const uint32_t *a , const uint32_t *b , int n ) {

		// Ensure this polynomial's capacity now to avoid reallocation
		ensureCapacity(n+1);

		// Initialize this polynomial by zero 'f(X) <- 0'
		setZero();

		SmallBinaryFieldPolynomial linear(getField());
		SmallBinaryFieldPolynomial l(getField());
		SmallBinaryFieldPolynomial li(getField());
		SmallBinaryFieldPolynomial tmp(getField());

		// Set the monomial 'linear(X) <- X'
		linear.setX();
		// Set the product 'l(X) <- (X-a[0])*(X-a[1])*...*(X-a[n-1])'
		l.buildFromRoots(a,n);
		// Ensure capacity now to avoid reallocation
		li.ensureCapacity(n+1);
		tmp.ensureCapacity(n+1);

		// Iterate over the Lagrange polynomials

		for ( int i = 0 ; i < n ; i++ ) {


			// ***************************************************************
			// *** BEGIN: Compute prefactor for 'i'-th Lagrange polynomial ***
			// ***************************************************************

			// 'u <- prod_{i!=j}(a[i]-a[j])
			uint32_t u = 1;
			for ( int j = 0 ; j < n ; j++ ) {
				if ( i != j ) {
					 u = this->gfPtr->mul(u,a[i]^a[j]);
				}
			}

			if ( u == 0 ) {
				cerr << "SmallBinaryFieldPolynomial::interpolate: "
				     << "Locators must be distinct." << endl;
				exit(EXIT_FAILURE);
			}

			// 'u <- (prod_{i!=j}(a[i]-a[j]))^(-1)'
			u = this->gfPtr->inv(u);

			// 'u <- b[i] * (\prod_{i!=j}(a[i]-a[j]))^(-1)'
			u = this->gfPtr->mul(u,b[i]);

			// ***************************************************************
			// **** END: 'u' is now the prefactor for the 'i'th Lagrange *****
			// **** polynomial ***********************************************
			// ***************************************************************


			// Determine the 'i'th Lagrange polynomial
			linear.setCoeff(0,this->gfPtr->neg(a[i]));
			divRem(li,tmp,l,linear);
			mul(li,li,u);

			// Determine the interpolation polynomial of the first 'i+1'
			// tuples '(a[0],b[0]),...,(a[i],b[i])' of degree '<= i'
			add(*this,*this,li);
		}
	}

	/**
	 * @brief
	 *            Multiplies the polynomial <i>g(X)</i> by the
	 *            monomial \f$X^n\f$ and stores the result
	 *            in <i>f</i>.
	 *
	 * @details
	 *            If <i>n>=0</i> the method lets
	 *            \f[
	 *             f(X)\leftarrow g(X)\cdot X^n.
	 *            \f]
	 *            If <i>n<0</i> the polynomial <i>f</i> will be set
	 *            as the quotient polynomial of <i>g</i> divided by
	 *            \f$X^n\f$ which is computed with the help of
	 *            <code>rightShift()</code>.
	 *
	 * @param f
	 *            Output polynomial.
	 *
	 * @param g
	 *            Input polynomial.
	 *
	 * @param n
	 *            The degree of the monomial that is multiplied with
	 *            <i>g</i>.
	 *
	 * @warning
	 *            If not sufficient memory could be provided, the method
	 *            prints an error message to <code>stderr</code> and exits
	 *            with status 'EXIT_FAILURE'.
	 */
	void SmallBinaryFieldPolynomial::leftShift
	( SmallBinaryFieldPolynomial & f ,
	  const SmallBinaryFieldPolynomial & g , int n ) {

		if ( n < 0 ) {
			rightShift(f,g,-n);
			return;
		}

		if ( !g.isZero() && n > 0 ) {

			int d = g.deg();

			f.ensureCapacity(d+n+1);
			f.degree = d+n;

			uint32_t *a , *b;
			a = f.coefficients + n+d;
			b = g.coefficients +  d;

			for ( int i = d ; i >= 0 ; i-- , a-- , b-- ) {
				*a = *b;
			}

			a = f.coefficients;
			for ( int i = 0 ; i < n ; i++ , a++ ) {
				*a = 0;
			}
		} else {
			f.assign(g);
		}
	}

	/**
	 * @brief
	 *            Computes the quotient of the polynomial <i>g(X)</i>
	 *            divided by the monomial \f$X^n\f$ and stores the result
	 *            in <i>f</i>.
	 *
	 * @details
	 *            see 'SmallBinaryFieldPolynomial.h'
	 */
	void SmallBinaryFieldPolynomial::rightShift
	( SmallBinaryFieldPolynomial & f ,
	  const SmallBinaryFieldPolynomial & g , int n ) {

		if ( n < 0 ) {
			leftShift(f,g,n);
			return;
		}

		f.assign(g);

		int m = f.deg();
		uint32_t *A = f.coefficients;

		for ( int i = 0 ; i+n <= m  ; i++ ) {
			A[i] = A[i+n];
		}

		if ( (f.degree -= n) < 0 ) {
			f.degree = -1L;
		}
	}

	/**
	 * @brief
	 *            Computes the bitwise exclusive or of the
	 *            <code>uint64_t</code> elements in
	 *            <code>a</code> with the elements in <code>b</code> and
	 *            stores the results in <code>c</code>
	 *
	 * @details
	 *            For \f$i=0,...,n-1\f$ the method performs
	 *            \f$c[i]\leftarrow a[i]\oplus b[i]</code> where \f$\oplus\f$
	 *            denotes the exclusive or operation on unsigned 64-bit
	 *            integers.
	 *
	 * @param c
	 *            output array that will contain the results of the xors of
	 *            the elements in <code>b</code> and the elements in
	 *            <code>a</code>
	 *
	 * @param a
	 *            first array of <code>uint64_t</code>'s
	 *
	 * @param b
	 *            second input array of <code>uint64_t</code>'s
	 *
	 * @param n
	 *            the number of processed <code>uint64_t</code>'s
	 *
	 * @warning
	 *            <code>a</code> or <code>b</code> do not contain at least
	 *            <code>n</code> valid <code>uint64_t</code>'s (i.e. if the
	 *            entries have not been set) or if <code>c</code>,
	 *            <code>a</code>, or <code>b</code> have not been allocated
	 *            to hold at least <code>n</code> <code>uint64_t</code>'s
	 *            calling this method may run into unexpected behavior.
	 */
	static void mxor64
	( register uint64_t *c ,
	  register uint64_t *a ,
	  register uint64_t *b , register int n ) {

		register int i;
		for ( i = 0 ; i < n ; i++ , c++ , a++ , b++ ) {
			*c = *a ^ *b;
		}
	}

	/**
	 * @brief
	 *            Computes the bitwise exclusive or of the
	 *            <code>uint32_t</code> elements in
	 *            <code>a</code> with the elements in <code>b</code> and
	 *            stores the results in <code>c</code>
	 *
	 * @details
	 *            For \f$i=0,...,n-1\f$ the method performs
	 *            \f$c[i]\leftarrow a[i]\oplus b[i]</code> where \f$\oplus\f$
	 *            denotes the exclusive or operation on unsigned 32-bit
	 *            integers.
	 *            <br><br>
	 *            This method uses the
	 *            \link mxor64(uint64_t*,uint64_t*,uint64_t*,int)\endlink
	 *            to achieve the desired effect for the first
	 *            \f$\lfloor n/2\rfloor\cdot 2\f$ elements.
	 *            If \f$n\equiv 1\mod 2\f$ the last entries are xored
	 *            separately. This is, in particular, beneficial for 64-bit
	 *            computers.
	 *
	 * @param c
	 *            output array that will contain the results of the xors of
	 *            the elements in <code>b</code> and the elements in
	 *            <code>a</code>
	 *
	 * @param a
	 *            first array of <code>uint32_t</code>'s
	 *
	 * @param b
	 *            second input array of <code>uint32_t</code>'s
	 *
	 * @param n
	 *            the number of processed <code>uint32_t</code>'s
	 *
	 * @warning
	 *            <code>a</code> or <code>b</code> do not contain at least
	 *            <code>n</code> valid <code>uint32_t</code>'s (i.e. if the
	 *            entries have not been set) or if <code>c</code>,
	 *            <code>a</code>, or <code>b</code> have not been allocated
	 *            to hold at least <code>n</code> <code>uint32_t</code>'s
	 *            calling this method may run into unexpected behavior.
	 */
	static inline void mxor32
	( uint32_t *c , uint32_t *a , uint32_t *b , int n ) {

		mxor64((uint64_t*)c,(uint64_t*)a,(uint64_t*)b,n/2);

		if ( n%2 ) {
			c[n-1] = a[n-1]^b[n-1];
		}
	}

	/**
	 * @brief
	 *            Computes the sum of two polynomials.
	 *
	 * @details
	 *            see 'SmallBinaryFieldPolynomial.h'
	 */
	void SmallBinaryFieldPolynomial::add
	( SmallBinaryFieldPolynomial & h ,
	  const SmallBinaryFieldPolynomial & f ,
	  const SmallBinaryFieldPolynomial & g ) {

		uint32_t *a , *b , *c;
		int m , n;

		// Obtain the higher degree 'n' and the lower
		// degree of the summands
		if ( f.deg() < g.deg() ) {
			m = f.degree;
			n = g.degree;
		} else {
			m = g.degree;
			n = f.degree;
		}

		// Ensure the output polynomial can store the sum
		h.ensureCapacity(n+1);

		// Ensure that 'a' contains the coefficients of the summand
		// that is of the higher degree 'n'; the coefficients of the
		// other polynomial of the lower degree 'm' are accessed by 'b'
		if ( f.deg() < g.deg() ) {
			b = f.coefficients;
			c = g.coefficients;
		} else {
			b = g.coefficients;
			c = f.coefficients;
		}

		// Temporarily, set the maximal possible degree
		h.degree = n;

		// 'a' accesses to coefficient of the output polynomial
		a = h.coefficients;

		// The first 'm+1' coefficients of the output polynomial
		// corresponds to the sum of the first respective 'm+1' elements
		// of the two input polynomials. Each of these sums, in turn,
		// correspond to a binary exclusive or operation
		mxor32(a,b,c,m+1);

		// The remaining coefficients simply are a copy of the corresponding
		// coefficient of the higher degree polynomial
		if ( a != c && n > m ) {
			memcpy(a+m+1,c+m+1,sizeof(uint32_t)*(n-m));
		}

		// The degree might be wrong, e.g., if the highest index coefficient
		// summed to zero. Therefore, ensure the degree is correct
		h.normalize();
	}

	/**
	 * @brief
	 *           Low-level method for multiplying two polynomials.
	 *
	 * @details
	 *           If
	 *           \f$f(X)=\sum_{i=0}^{m-1}a_iX^i\f$,
	 *           \f$g(X)=\sum_{j=0}^{n-1}b_jX^j\f$, and
	 *           \f$g(X)=\sum_{k=0}^{n+m-1}c_kX^k\f$ is the
	 *           product, then this method computes the \f$c_k\f$
	 *           from the <code>a[i]</code> and <code>b[j]</code>
	 *           and stores them in <code>c[k]</code>.
	 *
	 * @param c
	 *           output array which should be of size at least
	 *           <code>m+n-1</code>
	 *
	 * @param a
	 *           first input array which should contain at least
	 *           <code>m</code> valid elements in the specified
	 *           finite field
	 * @param m
	 *           number of valid elements in <code>a</code>
	 *
	 * @param b
	 *           second input array which should contain at least
	 *           <code>m</code> valid elements in the specified
	 *           finite field
	 *
	 * @param n
	 *           number of valid elements in <code>b</code>
	 *
	 * @param gf
	 *           pointer to the underlying binary finite field
	 *
	 * @warning
	 *           If <code>a</code> or <code>b</code> do not contain
	 *           valid representation of elements in the finite field
	 *           to where <code>gf</code> points or if <code>c</code>
	 *           was not allocated to store at least <code>m+n-1</code>
	 *           elements the method may run into unexpected behavior.
	 *           Furthermore, the entries of <code>a</code>,
	 *           <code>b</code>, and <code>c</code> should not cross.
	 *
	 */
	static void TradMul
	( register uint32_t *c ,
	  register uint32_t *a , register int m ,
	  register uint32_t *b , register int n ,
	  const SmallBinaryField *gfPtr ) {

		register uint32_t *d , *e;
		register int l = m+n-1;

		d = c;

		if ( l > 0 )
			memset(d,0,sizeof(uint32_t)*l);

        register uint32_t carry = 0;

		for ( register int i = 0 ; i < m ; i++ , c++ , a++ ) {

			d = c;
			e = b;

			for ( int j = 0 ; j < n ; j++ , d++ , e++ ) {

                carry = gfPtr->mul(*a,*e);

				*d ^= carry;
			}
		}
	}

	/**
	 * @brief
	 *            Computes the product of two polynomials.
	 *
	 * @details
	 *            see 'SmallBinaryFieldPolynomial.h'
	 *
	 */
	void SmallBinaryFieldPolynomial::mulUncheck
	( SmallBinaryFieldPolynomial & h ,
	  const SmallBinaryFieldPolynomial & f ,
	  const SmallBinaryFieldPolynomial & g ) {

		uint32_t *a , *b , *c;
		int m , n;

		if ( f.deg() < g.deg() ) {
			m = f.deg();
			n = g.deg();
			a = f.coefficients;
			b = g.coefficients;
		} else {
			m = g.deg();
			n = f.deg();
			a = g.coefficients;
			b = f.coefficients;
		}

        h.ensureCapacity(m+n+1);
		h.degree = m+n;

		c = h.coefficients;

		TradMul(c,a,m+1,b,n+1,h.gfPtr);

		h.normalize();
	}

	/**
	 * @brief
	 *            Computes the product of two polynomials.
	 *
	 * @details
	 *            see 'SmallBinaryFieldPolynomial.h'
	 */
	void SmallBinaryFieldPolynomial::mul
	( SmallBinaryFieldPolynomial & h ,
	  const SmallBinaryFieldPolynomial & f ,
	  const SmallBinaryFieldPolynomial & g ) {

		if ( h.gfPtr->getDefiningPolynomial().rep != f.gfPtr->getDefiningPolynomial().rep ||
			 h.gfPtr->getDefiningPolynomial().rep != g.gfPtr->getDefiningPolynomial().rep ) {
			cerr << "SmallBinaryFieldPolynomial::mul: Polynomial's "
					"corresponding finite fields must be of same reference."
				 << endl;
			exit(EXIT_FAILURE);
		}

		if ( f.isZero() || g.isZero() ) {
			h.setZero();
			return;
		}

		if ( &h == &f && &h == &g ) {
			SmallBinaryFieldPolynomial tf(f) , tg(g);
			mulUncheck(h,tf,tg);
		} else if ( &h == &f ) {
			SmallBinaryFieldPolynomial tf(f);
			mulUncheck(h,tf,g);
		} else if ( &h == &g ) {
			SmallBinaryFieldPolynomial tg(g);
			mulUncheck(h,f,tg);
		} else {
			mulUncheck(h,f,g);
		}
	}

	/**
	 * @brief
	 *           Low-level method for the multiplication of a polynomial
	 *           by a scalar.
	 *
	 * @details
	 *           This method multiplies the first <code>n</code> entries of
	 *           <code>a</code> by <code>s</code> inplace. More precisely,
	 *           \f$a[i]\leftarrow a[i]\cdot s\f$ for \f$i=0,...,n-1\f$.
	 *           The multiplication is performed in the field to where
	 *           <code>gfPtr</code> points.
	 *
	 * @param a
	 *           Field of <code>n</code> valid elements in the field pointed
	 *           by <code>gfPtr</code>.
	 *
	 * @param s
	 *           Scalar which is an element of the field to where
	 *           <code>gfPtr</code> points.
	 *
	 * @param gfPtr
	 *           Pointer to an initialized small binary field.
	 *
	 * @warning
	 *           If <code>a</code> contains invalid elements in the field
	 *           to where <code>gfPtr</code> points, or if <code>s</code>
	 *           is an element field element, or if <code>gfPtr</code>
	 *           does not point to a successfully initialized small binary
	 *           field, the method may run into undocumented behavior.
	 */
	static void mul
	( register uint32_t *a , register uint32_t s , register int n ,
	  const SmallBinaryField *gfPtr ) {

		register int i;

		for ( i = 0 ; i < n ; i++ , a++ ) {
			*a = gfPtr->mul(*a,s);
		}

	}

	/**
	 * @brief
	 *            Multiplies a polynomial with a constant element of the
	 *            finite field.
	 *
	 * @details
	 *            see 'SmallBinaryFieldPolynomial.h'
	 */
	void SmallBinaryFieldPolynomial::mul
	( SmallBinaryFieldPolynomial & h ,
	  const SmallBinaryFieldPolynomial & f ,
	  uint32_t s ) {

		// Copy 'f' to 'h'
		h = f;

		// Multiply the coefficients inplace
		thimble::mul(h.coefficients,s,h.degree+1,h.gfPtr);
	}

	/**
	 * @brief
	 *            Performs Euclidean division of two polynomial to obtain
	 *            its quotient and remainder.
	 *
	 * @details
	 *            see 'SmallBinaryFieldPolynomial.h'
	 */
	void SmallBinaryFieldPolynomial::divRem
	( SmallBinaryFieldPolynomial & q ,
	  SmallBinaryFieldPolynomial & r ,
	  const SmallBinaryFieldPolynomial & a ,
	  const SmallBinaryFieldPolynomial & b ) {

		const SmallBinaryField *gfPtr = q.gfPtr;

		if ( r.gfPtr->getDefiningPolynomial().rep != gfPtr->getDefiningPolynomial().rep ||
			 a.gfPtr->getDefiningPolynomial().rep != gfPtr->getDefiningPolynomial().rep ||
			 b.gfPtr->getDefiningPolynomial().rep != gfPtr->getDefiningPolynomial().rep ) {
			cerr << "SmallBinaryFieldPolynomial::divRem: "
				 << "The given polynomial's underlying finite field must "
				 << "be of the same reference." << endl;
			exit(EXIT_FAILURE);
		}

		if ( &q == &r ) {
			cerr << "SmallBinaryFieldPolynomial::divRem: "
			     << "Quotient and remainder polynomial must be of different "
			     << "references." << endl;
			exit(EXIT_FAILURE);
		}

		if ( &q == &b  || &r == &b ) {
			SmallBinaryFieldPolynomial tb(b);
			divRem(q,r,a,tb);
			return;
		}

		if ( b.isZero() ) {
			cerr << "SmallBinaryFieldPolynomial::divRem: Division by zero."
			<< endl;
			exit(EXIT_FAILURE);
		}

		int m , n;

		m = b.deg();
		n = a.deg();

		uint32_t u;

		r = a;
		u = b.coefficients[m];
		u = gfPtr->inv(u);

		int d = n-m;
		if ( d < 0 ) {
			d = -1;
		}

		q.setZero();
		q.ensureCapacity(d+1);
		q.degree = d;

		uint32_t *Q , *R , *B;

		Q = q.coefficients+d;
		for ( int i = n-m ; i >= 0 ; i-- , Q-- ) {

			if ( r.deg() == m+i ) {

				R = r.coefficients+i;
				B = b.coefficients;

				*Q = gfPtr->mul(R[m],u);

				for ( int j = 0 ; j <= m ; j++ , B++ , R++ ) {
					*R ^= gfPtr->mul(*Q,*B);
				}

				r.normalize();

			} else {
				*Q = 0;
			}
		}
	}

    /**
     * @brief
     *          Computes the composition of two polynomials with
     *          coefficients in a small binary field.
     *
     * @details
     *          Given two polynomials <i>f</i> and <i>g</i> the composition
     *          <i>h=f(g)</i> is computed.
     *
     * @param h
     *          The composition of the two polynomials <i>f</i> and <i>g</i>
     *          on input, i.e., <i>f(g)</i>.
     *
     * @param f
     *          The outer polynomial.
     *
     * @param g
     *          The inner polynomial.
     *
     * @warning
     *          If not enough memory could be provided, an error
     *          message is printed to <code>stderr</code> and the
     *          program exits with status 'EXIT_FAILURE'.
     */
    void SmallBinaryFieldPolynomial::eval
    ( SmallBinaryFieldPolynomial & h ,
      const SmallBinaryFieldPolynomial & f ,
      const SmallBinaryFieldPolynomial & g ) {

        if ( &h == &f ) {
            SmallBinaryFieldPolynomial tf(f);
            eval(h,tf,g);
            return;
        }

        if ( &h == &g ) {
            SmallBinaryFieldPolynomial tg(g);
            eval(h,f,tg);
            return;
        }

        SmallBinaryFieldPolynomial tmp(f.getField());
        tmp.ensureCapacity((f.degree*(g.degree-1)));

        h.setZero();
        int d = f.deg();

        // Run Horner's method
        for ( int j = d ; j >= 0 ; j--  ) {
            mul(tmp,h,g);
            swap(h,tmp);
            h.setCoeff(0,f.getCoeff(j)^h.getCoeff(0));
        }
    }

	/**
	 * @brief
	 *            Computes partial greatest common divisor of two
	 *            polynomials.
	 *
	 * @details
	 *            see 'SmallBinaryFieldPolynomial.h'
	 */
	void SmallBinaryFieldPolynomial::pgcd
	( SmallBinaryFieldPolynomial & g ,
	  SmallBinaryFieldPolynomial & s ,
	  SmallBinaryFieldPolynomial & t ,
	  const SmallBinaryFieldPolynomial & a ,
	  const SmallBinaryFieldPolynomial & b ,
	  int d ) {

		SmallBinaryFieldPolynomial
			q(a.getField()) ,
			r0(a.getField()) ,
			s0(a.getField()) ,
			t0(a.getField()) ,
			r1(a.getField()) ,
			s1(a.getField()) ,
			t1(a.getField()) ,
			tmp(a.getField());

	    r0 = a; s0.setOne() ; t0.setZero();
	    r1 = b; s1.setZero(); t1.setOne();

	    while ( !r1.isZero() && r0.deg() >= d ) {

	    	divRem(q,tmp,r0,r1);

	    	r0.swap(r1);
	    	s0.swap(s1);
	    	t0.swap(t1);

	    	mul(tmp,q,r0);
	    	sub(r1,r1,tmp);

	    	mul(tmp,q,s0);
	    	sub(s1,s1,tmp);

	    	mul(tmp,q,t0);
	    	sub(t1,t1,tmp);
	    }

	    s = s0;
	    t = t0;
	    g = r0;
	}

	/**
	 * @brief
	 *            Applies the extended Euclidean algorithm to
	 *            two polynomials and stores the sequence of
	 *            the algorithm in vectors.
	 *
	 * @details
	 *            The method initializes
	 *            <pre>
	 *             g[0] = s[0]*a+t[0]*b
	 *             g[1] = s[1]*a+t[1]*b
	 *            </pre>
	 *            where <code>s[0]=1</code>, <code>t[0]=0</code>
	 *            and <code>s[1]=0</code>, and <code>t[1]=1</code>.
	 *            Then for each <code>i>=2</code> the method computes
	 *            <pre>
	 *             q[i] = g[i-2] rem g[i-1]
	 *            </pre>
	 *            and stores
	 *            <pre>
	 *             g[i] = g[i-2] - q[i]*g[i-1]
	 *             s[i] = s[i-2] - q[i]*s[s-1]
	 *             t[i] = t[i-2] - q[i]*t[i-1]
	 *            </pre>
	 *            in the output vectors passed to this method until
	 *            <code>g.back()</code> is a greatest common divisor
	 *            of <code>a</code> and <code>b</code>. In this way,
	 *            the vectors <code>g, s</code> and <code>t</code>
	 *            have the same length <code>n</code> and for each
	 *            <code>i=0,...,n-1</code> the relation
	 *            <pre>
	 *             g[i] = s[i]*a+t[i]*b
	 *            </pre>
	 *            holds.
	 *
	 * @attention
	 *            Note that, if <code>a</code> and <code>b</code> are
	 *            co-prime, the value <code>g.back()</code> does not
	 *            necessarily need to be equals 1 but is a non-zero
	 *            constant, i.e., <code>g.back.deg()==0</code>.
	 *
	 * @param g
	 *            Output vector of polynomials; see details.
	 *
	 * @param s
	 *            Output vector of polynomials; see details.
	 *
	 * @param t
	 *            Output vector of polynomials; see details.
	 *
	 * @param a
	 *            First input polynomial.
	 *
	 * @param b
	 *            Second input polynomial.
	 *
	 * @warning
	 *            If not enough memory could be provided, an
	 *            error message is printed to <code>stderr</code>
	 *            and the program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If <code>g</code>, <code>s</code> or <code>t</code>
	 *            are vectors containing <code>a</code> or <code>b</code>
	 *            or if <code>g</code>, <code>s</code> and <code>t</code>
	 *            are not of pair-wise different reference, this method
	 *            runs into undocumented behavior.
	 */
	void SmallBinaryFieldPolynomial::xgcd
	( std::vector<SmallBinaryFieldPolynomial> & r ,
	  std::vector<SmallBinaryFieldPolynomial> & s ,
	  std::vector<SmallBinaryFieldPolynomial> & t ,
	  const SmallBinaryFieldPolynomial & a ,
	  const SmallBinaryFieldPolynomial & b ) {


		SmallBinaryFieldPolynomial
			q(a.getField()) ,
			r0(a.getField()) ,
			s0(a.getField()) ,
			t0(a.getField()) ,
			r1(a.getField()) ,
			s1(a.getField()) ,
			t1(a.getField()) ,
			tmp(a.getField());

		r.clear();
		s.clear();
		t.clear();

	    r0 = a; s0.setOne() ; t0.setZero();
	    r1 = b; s1.setZero(); t1.setOne();

	    r.push_back(r0);
	    s.push_back(s0);
	    t.push_back(t0);

	    while ( !r1.isZero() && r0.deg() >= 0 ) {

	    	r.push_back(r1);
	    	s.push_back(s1);
	    	t.push_back(t1);

	    	divRem(q,tmp,r0,r1);

	    	r0.swap(r1);
	    	s0.swap(s1);
	    	t0.swap(t1);

	    	mul(tmp,q,r0);
	    	sub(r1,r1,tmp);

	    	mul(tmp,q,s0);
	    	sub(s1,s1,tmp);

	    	mul(tmp,q,t0);
	    	sub(t1,t1,tmp);
	    }
	}

	/**
	 * @brief
	 *            Applies the extended Euclidean algorithm to two
	 *            polynomials.
	 *
	 * @details
	 *            This method computes polynomials <code>g</code>,
	 *            <code>s</code> and <code>t</code> such that
	 *            <code>g</code> is a greatest common divisor of
	 *            <code>a</code> and <code>b</code> where
	 *            <pre>
	 *             g = s * a + t *b.
	 *            </pre>
	 *
	 * @attention
	 *            Note that, if <code>a</code> and <code>b</code> are
	 *            co-prime, the value <code>g</code> does not
	 *            necessarily need to be equals 1 but is a non-zero
	 *            constant, i.e., <code>g.deg()==0</code>.
	 *
	 * @param g
	 *            Output polynomial that is a greatest common divisor
	 *            of <code>a</code> and <code>b</code>.
	 *
	 * @param s
	 *            Output polynomial; see details above.
	 *
	 * @param t
	 *            Output polynomial; see details above.
	 *
	 * @param a
	 *            First input polynomial.
	 *
	 * @param b
	 *            Second input polynomial.
	 *
	 * @warning
	 *            If not enough memory could be provided, an
	 *            error message is printed to <code>stderr</code>
	 *            and the program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If <code>g</code>, <code>s</code> and <code>t</code>
	 *            are not of pair-wise different reference, this method
	 *            runs into undocumented behavior.
	 */
	void SmallBinaryFieldPolynomial::xgcd
	( SmallBinaryFieldPolynomial & g ,
	  SmallBinaryFieldPolynomial & s ,
	  SmallBinaryFieldPolynomial & t ,
	  const SmallBinaryFieldPolynomial & a ,
	  const SmallBinaryFieldPolynomial & b ) {

		SmallBinaryFieldPolynomial
			q(a.getField()) ,
			r0(a.getField()) ,
			s0(a.getField()) ,
			t0(a.getField()) ,
			r1(a.getField()) ,
			s1(a.getField()) ,
			t1(a.getField()) ,
			tmp(a.getField());

	    r0 = a; s0.setOne() ; t0.setZero();
	    r1 = b; s1.setZero(); t1.setOne();

	    while ( !r1.isZero() ) {

	    	divRem(q,tmp,r0,r1);

	    	r0.swap(r1);
	    	s0.swap(s1);
	    	t0.swap(t1);

	    	mul(tmp,q,r0);
	    	sub(r1,r1,tmp);

	    	mul(tmp,q,s0);
	    	sub(s1,s1,tmp);

	    	mul(tmp,q,t0);
	    	sub(t1,t1,tmp);
	    }

	    s = s0;
	    t = t0;
	    g = r0;
	}

	/**
	 * @brief
	 *            Computes a greatest common divisor of of two
	 *            polynomials.
	 *
	 * @attention
	 *            Note that, if <code>a</code> and <code>b</code> are
	 *            co-prime, the value <code>g</code> does not
	 *            necessarily need to be equals 1 but is a non-zero
	 *            constant, i.e., <code>g.deg()==0</code>.
	 *
	 * @param g
	 *            Output polynomial that represents a greatest common
	 *            divisor of <code>a</code> and <code>b</code>.
	 *
	 * @param a
	 *            First input polynomial.
	 *
	 * @param b
	 *            Second input polynomial.
	 *
	 * @warning
	 *            If not enough memory could be provided, an
	 *            error message is printed to <code>stderr</code>
	 *            and the program exits with status 'EXIT_FAILURE'.
	 */
	void SmallBinaryFieldPolynomial::gcd
	( SmallBinaryFieldPolynomial & g ,
	  const SmallBinaryFieldPolynomial & a ,
	  const SmallBinaryFieldPolynomial & b ) {

		SmallBinaryFieldPolynomial s(a.getField()) , t(a.getField());

		SmallBinaryFieldPolynomial::xgcd(g,s,t,a,b);
	}

	/**
	 * @brief
	 *            Checks whether two polynomials are co-prime.
	 *
	 * @details
	 *            Two polynomials <code>a</code> and <code>b</code>
	 *            are said to be co-prime if the only non-zero
	 *            polynomials <code>c</code> that divide <code>a</code>
	 *            and <code>b</code> are constants.
	 *
	 * @param a
	 *            First polynomial.
	 *
	 * @param b
	 *            Second polynomial.
	 *
	 * @return
	 *            <code>true</code> if <code>a</code> and <code>b</code>
	 *            are co-prime and, otherwise, <code>false</code>.
	 *
	 * @warning
	 *            If not enough memory could be provided, an
	 *            error message is printed to <code>stderr</code>
	 *            and the program exits with status 'EXIT_FAILURE'.
	 */
	bool SmallBinaryFieldPolynomial::areCoprime
	( const SmallBinaryFieldPolynomial & a , const SmallBinaryFieldPolynomial & b ) {

		SmallBinaryFieldPolynomial g(a.getField());

		SmallBinaryFieldPolynomial::gcd(g,a,b);

		return g.deg()==0;
	}

	/**
	 * @brief
	 *            Computes the modular inverse of a polynomial.
	 *
	 * @details
	 *            Given two polynomials \f$x\f$ and \f$m\f$, this function
	 *            determines a polynomial \f$y\f$ with
	 *            \f$\deg(y)<\deg(m)\f$ such that
	 *            \f$(y\cdot x)~mod~m=1\f$. If such a polynomial exists,
	 *            this method will succeed in finding such a polynomial
	 *            and returns <code>true</code>; otherwise, if such a
	 *            polynomial does not exist, <code>y</code> will be
	 *            left unchanged and the function returns
	 *            <code>false</code>.
	 *
	 * @param y
	 *            Output polynomial that forms the inverse of <code>x</code>
	 *            modulo <code>m</code>.
	 *
	 * @param x
	 *            Input polynomial of which the modular inverse is
	 *            computed.
	 *
	 * @param m
	 *            Input modulus polynomial.
	 *
	 * @warning
	 *            If not enough memory could be provided, an
	 *            error message is printed to <code>stderr</code>
	 *            and the program exits with status 'EXIT_FAILURE'.
	 */
	bool SmallBinaryFieldPolynomial::invMod
	( SmallBinaryFieldPolynomial & y , const SmallBinaryFieldPolynomial & x ,
	  const SmallBinaryFieldPolynomial & m ) {

		if ( &y == &m ) {
			return invMod(y,x,SmallBinaryFieldPolynomial(m));
		}

		SmallBinaryFieldPolynomial g(x.getField()) , s(x.getField()) , t(x.getField());

		SmallBinaryFieldPolynomial::xgcd(g,s,t,x,m);

		if ( g.deg() != 0 ) {
			return false;
		}

		mul(y,s,x.getField().inv(g.getCoeff(0)));

		rem(y,y,m);

		return true;
	}

	/**
	 * @brief
	 *           Prints a text representation of a polynomial to the specified
	 *           output stream.
	 *
	 * @details
	 *           see 'SmallBinaryFieldPolynomial.h'
	 */
	ostream & operator<<
			( ostream & out , const SmallBinaryFieldPolynomial & f ) {

	    int d = f.deg();

	    out << std::uppercase;

	    if ( d >= 0 ) {

	    	// If the degree is larger than 0, then at least
	    	// one element is non-zero causing a non-empty printout

	    	for ( int i = 0 ; i <= d ; i++ ) {

	    		if ( f.getCoeff(i) != 0) {

	    			out << "0x" << std::hex << f.getCoeff(i);
	    			if ( i > 0 ) {
	    				out << "*X";
	    				if ( i > 1 )
	    					out << "^" << std::dec << i;
	    			}

	    			if ( i < d )
	    				out << "+";
	    		}
	      }

	    } else {

	    	// If the polynomial is zero, print a hexadecimal zero
	    	out << "0x0";
	    }

	    out << std::dec;

		return out;
	}
}



