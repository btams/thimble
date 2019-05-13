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
 * @file SmallBinaryPolynomial.cpp
 *
 * @brief
 *            Implements the functionalities that are provided by
 *            'SmallBinaryPolynomial.h' which are related with a class to
 *            represent and compute with polynomials of small degree having
 *            boolean coefficients.
 *
 * @author Benjamin Tams
 *
 * @see thimble::SmallBinaryPolynomial
 */

#include <cstdlib>
#include <iostream>

#include <thimble/math/numbertheory/SmallBinaryPolynomial.h>

using namespace std;

/**
 * @brief The library's namespace
 */
namespace thimble {

	/**
	 * @brief
	 *            Performs Euclidean division of two polynomial to obtain
	 *            its quotient and remainder.
	 *
	 * @details
	 *            The Euclidean quotient \f$q(X)\f$ and the
	 *            Euclidean remainder \f$r(X)\f$ of two polynomials
	 *            \f$a(X)\f$ and \f$b(X)\f$ are defined to fulfill the
	 *            relation \f$q(X)\cdot b(X)=a(X)+r(X)\f$ where
	 *            \f$\deg(r(X))<\deg(b(X))\f$. If \f$b(X)\neq 0\f$ then
	 *            the Euclidean quotient and remainder are known to exist
	 *            and that they are unique.
	 *
	 * @warning
	 *            If <code>b</code> is the zero polynomial or if
	 *            <code>q</code> and <code>r</code> are of the same
	 *            reference then an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *             -1.
	 *
	 * @param a
	 *            numerator
	 *
	 * @param b
	 *            denominator
	 *
	 * @param q
	 *            will contain the quotient of the Euclidean division of
	 *            <code>a</code> divided by <code>b</code>
	 *
	 * @param r
	 *            will contain the remainder of the Euclidean division
	 *            of <code>a</code> divided by <code>b</code>
	 */
	void SmallBinaryPolynomial::divRem
	( SmallBinaryPolynomial & q , SmallBinaryPolynomial & r ,
	  const SmallBinaryPolynomial & a , const SmallBinaryPolynomial & b ) {

		if ( &q == &r ) {
			cerr << "SmallBinaryPolynomial::divRem: Quotient and remainder "
				 << " must be of different references." << endl;
			exit(EXIT_FAILURE);
		}

		if ( b.rep == 0 ) {
			cerr << "SmallBinaryPolynomial::divRem: Division by zero" << endl;
			exit(EXIT_FAILURE);
		}

		SmallBinaryPolynomial c = b;

		int n , m;
		n = a.deg();
		m = b.deg();

		r.rep = a.rep;
		q.rep = 0;

		for ( int i = n-m ; i >= 0 ; i--) {

			if ( r.deg() == m+i ) {
				q.setCoeff(i);
				r.rep ^= c.rep<<i;
			}

		}
	}

	/**
	 * @brief
	 *            Performs the Extended Euclidean Algorithm.
	 *
	 * @details
	 *            Computes <code>g</code>, <code>s<code>, and
	 *            <code>t</code> such that
	 *            <code>g=s*a+t*b</code> where <code>g</code>
	 *            is the greatest common divisor of
	 *            <code>a</code> and <code>b</code>.
	 *            <br><br>
	 *            The greatest common divisor <code>g</code>of two
	 *            polynomials <code>a</code> and <code>b</code> is
	 *            defined to be the polynomial of highest degree
	 *            with <code>g</code> divides <code>a</code> and
	 *            <code>g</code> divides <code>b</code>. For general
	 *            polynomials the greatest divisor is unique up
	 *            to a constant factor. In case of binary polynomials
	 *            it is even unique.
	 *
	 * @param g
	 *            Will contain the greatest common divisor of
	 *            <code>a</code> and <code>b</code>
	 *
	 * @param s
	 *            see details
	 *
	 * @param t
	 *            see details
	 *
	 * @param a
	 *            first argument
	 *
	 * @param b
	 *            second argument
	 */
	void SmallBinaryPolynomial::xgcd
	( SmallBinaryPolynomial & g ,
	  SmallBinaryPolynomial & s , SmallBinaryPolynomial & t ,
	  const SmallBinaryPolynomial & a , const SmallBinaryPolynomial & b) {

		SmallBinaryPolynomial r0 , s0 , t0 , r1 , s1 , t1 , r2 , s2 , t2 , q;

		r0 = a; s0 = 1; t0 = 0;
		r1 = b; s1 = 0; t1 = 1;

		while ( r1.rep ) {

			divRem(q,r2,r0,r1);
			s2 = s0.sub(q.mul(s1));
			t2 = t0.sub(q.mul(t1));

			r0 = r1;
			s0 = s1;
			t0 = t1;

			r1 = r2;
			s1 = s2;
			t1 = t2;
		}

		g = r0;
		s = s0;
		t = t0;
	}

	/**
	 * @brief Selects an irreducible polynomial of specified degree
	 *        from a precomputed list.
	 *
	 * @details
	 *        The degree must be between 0 and 30. Otherwise an error
	 *        message is printed to <code>stderr</code> and the program
	 *        exits with status 'EXIT_FAILURE'.
	 *        <br><br>
	 *        The result will be a polynomial that is irreducible and
	 *        that is of the specified degree. The function simply selects
	 *        the polynomial from a fixed list in where for each degree
	 *        (between 1 and 30) there is exactly one polynomial. This
	 *        is in particular useful for computing with finite fields.
	 *
	 * @warning
	 *        If degree is not between 1 and 31 an error message is
	 *        printed to <code>stderr</code> and the program exits with
	 *        status 'EXIT_FAILURE'.
	 *
	 * @param degree
	 *             The degree of the sought polynomial
	 *
	 * @return     An irreducible polynomial of specified degree
	 */
	SmallBinaryPolynomial SmallBinaryPolynomial::irreducible( int degree ) {

		if ( degree < 1 || degree > 30 ) {
			cerr << "The degree of the irreducible polynomial must range between "
				 << "1 and 30" << endl;
			exit(EXIT_FAILURE);
		}

		switch(degree) {
		case 1:
			return SmallBinaryPolynomial(2);
		case 2:
			return SmallBinaryPolynomial(7);
		case 3:
			return SmallBinaryPolynomial(11);
		case 4:
			return SmallBinaryPolynomial(19);
		case 5:
			return SmallBinaryPolynomial(37);
		case 6:
			return SmallBinaryPolynomial(67);
		case 7:
			return SmallBinaryPolynomial(131);
		case 8:
			return SmallBinaryPolynomial(283);
		case 9:
			return SmallBinaryPolynomial(515);
		case 10:
			return SmallBinaryPolynomial(1033);
		case 11:
			return SmallBinaryPolynomial(2053);
		case 12:
			return SmallBinaryPolynomial(4105);
		case 13:
			return SmallBinaryPolynomial(8219);
		case 14:
			return SmallBinaryPolynomial(16417);
		case 15:
			return SmallBinaryPolynomial(32771);
		case 16:
			return SmallBinaryPolynomial(65579);
		case 17:
			return SmallBinaryPolynomial(131081);
		case 18:
			return SmallBinaryPolynomial(262153);
		case 19:
			return SmallBinaryPolynomial(524327);
		case 20:
			return SmallBinaryPolynomial(1048585);
		case 21:
			return SmallBinaryPolynomial(2097157);
		case 22:
			return SmallBinaryPolynomial(4194307);
		case 23:
			return SmallBinaryPolynomial(8388641);
		case 24:
			return SmallBinaryPolynomial(16777243);
		case 25:
			return SmallBinaryPolynomial(33554441);
		case 26:
			return SmallBinaryPolynomial(67108891);
		case 27:
			return SmallBinaryPolynomial(134217767);
		case 28:
			return SmallBinaryPolynomial(268435459);
		case 29:
			return SmallBinaryPolynomial(536870917);
		case 30:
			return SmallBinaryPolynomial(1073741827);
		default:
			break;
		}

		// For the compiler even though this is not reachable
		return SmallBinaryPolynomial();
	}
}



