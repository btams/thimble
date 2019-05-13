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
 * @file SmallBinaryPolynomial.h
 *
 * @brief
 *            Provides a class to represent and compute with polynomials
 *            of small degree having boolean coefficients.
 *
 * @author Benjamin Tams
 *
 * @see thimble::SmallBinaryPolynomial
 */

#ifndef THIMBLE_SMALLBINARYPOLYNOMIAL_H_
#define THIMBLE_SMALLBINARYPOLYNOMIAL_H_

#include <stdint.h>
#include <iostream>

#include <thimble/dllcompat.h>
#include <thimble/math/MathTools.h>

/**
 * @brief The library's namespace
 */
namespace thimble {

	/**
	 * @brief Represents polynomials of small degree having binary
	 *        coefficients, i.e. either 1 (<code>true</code>) or
	 *        0 (<code>false</code>).
	 */
	class THIMBLE_DLL SmallBinaryPolynomial {

	public:
		/**
		 * @brief Representation of the polynomial.
		 *
		 * @details
		 *        The <i>i</i>-th order bit of \link rep\endlink corresponds to
		 *        the <i>i</i>-th coefficient of the polynomial. More precisely,
		 *        the integer \f$rep=\sum_{i=0}^db_i2^i\f$ represents the
		 *        polynomial \f$f(X)=\sum_{i=0}^db_iX^i\f$ in the indeterminate
		 *        <i>X</i> and vice versa. If furthermore <i>d</i> is such that
		 *        \f$ b_i\neq 0\f$ then <i>d</i> is equals the degree of the
		 *        polynomial.
		 */
		uint64_t rep;

		/**
		 * @brief   Creates a binary polynomial represented bit-wisely by the
		 *          specified integer
		 *
		 * @param rep
		 *          integer representing the created polynomial
		 *
		 */
		inline SmallBinaryPolynomial( uint64_t rep = 0 ) {
			this->rep = rep;
		}

		/**
		 * @brief Copy constructor
		 *
		 * @details
		 *             Creates an instance of the
		 *             \link SmallBinaryPolynomial\endlink class that
		 *             is a copy of the given polynomial.
		 *
		 * @param f
		 *             the polynomial of which the copy is created
		 */
		inline SmallBinaryPolynomial( const SmallBinaryPolynomial & f ) {
			this->rep = f.rep;
		}

		/**
		 * @brief
		 *         Assignment operator
		 *
		 * @details
		 *         Sets this polynomial to a copy of the given polynomial.
		 *
		 * @param f
		 *         The polynomial of which this polynomial will become a
		 *         copy of.
		 */
		inline SmallBinaryPolynomial & operator=
		( const SmallBinaryPolynomial & f ) {
			this->rep = f.rep; return *this;
		}

		/**
		 * @brief     Accesses the degree of the polynomial
		 *
		 * @details
		 *            The degree of a polynomial \f$f(X)\f$ is the least
		 *            integer <i>d</i> such that \f$f(X)=\sum_{i=0}^df_iX^i\f$
		 *            with \f$f_d\neq 0\f$. In case, the polynomial
		 *            is identical to zero the degree is -1
		 *
		 * @return
		 *            The degree of the polynomial
		 */
		inline int deg() const {
			return MathTools::numBits(this->rep)-1;
		}

		/**
		 * @brief     Access the <i>i</i>-th coefficient of the polynomial
		 *
		 * @details
		 *            The <i>i</i>-th coeffient of a polynomial \f$f(X)\f$
		 *            is the value \f$f_i\in\{0,1\}\f$ where
		 *            \f$f(X)=\sum_{j\geq 0}f_jX^i\f$.
		 *            <br><br>
		 *            Note that the result will be <code>0</code> if
		 *            <code>i<0</code> or if <code>i>63</code>.
		 *
		 * @param i
		 *            The index of the coefficient that is accessed
		 *
		 * @return
		 *            The <i>i</i>-th coefficient of the polynomial
		 *
		 * @see setCoeff(int)
		 * @see clearCoeff(int)
		 * @see setCoeff(int,bool)
		 */
		inline bool getCoeff( int i ) const {
			if ( i < 0 ) {
				return false;
			}
			return (this->rep&(((uint64_t)1)<<i))!=0?true:false;
		}

		/**
		 * @brief     Sets the <i>i</i>-th coefficient of the polynomial as
		 *            1 (true).
		 *
		 * @details
		 *            The <i>i</i>-th coeffient of a polynomial \f$f(X)\f$
		 *            is the value \f$f_i\in\{0,1\}\f$ where
		 *            \f$f(X)=\sum_{j\geq 0}f_jX^i\f$.
		 *
		 * @warning
		 *            If <code>i<0</code> or if <code>i>63</code> an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 *
		 * @param i
		 *            The index of the accessed coefficient
		 *
		 * @see getCoeff(int)
		 * @see clearCoeff(int)
		 * @see setCoeff(int,bool)
		 */
		inline void setCoeff( int i ) {
			if ( i < 0 || i > 63 ) {
				std::cerr <<
				"Accessed coefficient must range between 0 and 63"
				<< std::endl;
				exit(EXIT_FAILURE);
			}
			this->rep |= (((uint64_t)1)<<i);
		}

		/**
		 * @brief     Unsets the <i>i</i>-th coefficient of the polynomial as
		 *            0 (false).
		 *
		 * @details
		 *            The <i>i</i>-th coeffient of a polynomial \f$f(X)\f$
		 *            is the value \f$f_i\in\{0,1\}\f$ where
		 *            \f$f(X)=\sum_{j\geq 0}f_jX^i\f$.
		 *
		 * @warning
		 *            If <code>i<0</code> or if <code>i>63</code> an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 *
		 * @param i
		 *            The index of the accessed coefficient
		 *
		 * @see getCoeff(int)
		 * @see setCoeff(int)
		 * @see setCoeff(int,bool)
		 */
		inline void clearCoeff( int i ) {
			if ( i < 0 || i > 63 ) {
				std::cerr <<
				"Accessed coefficient must range between 0 and 63"
				<< std::endl;
				exit(EXIT_FAILURE);
			}
			this->rep &= ~(((uint64_t)1)<<i);
		}

		/**
		 * @brief     Sets the <i>i</i>-th coefficient of the polynomial as
		 *            <code>c</code> which is either 1 (true) or 0 (false).
		 *
		 * @details
		 *            The <i>i</i>-th coefficient of a polynomial \f$f(X)\f$
		 *            is the value \f$f_i\in\{0,1\}\f$ where
		 *            \f$f(X)=\sum_{j\geq 0}f_jX^i\f$.
		 *
		 * @warning
		 *            If <code>i<0</code> or if <code>i>63</code> an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 *
		 * @param i
		 *            The index of the accessed coefficient
		 *
		 * @param c
		 *            The new value of the accessed coefficient which is
		 *            either 1 (<code>true</code>) or 0 (<code>false</code>)
		 *
		 * @see getCoeff(int)
		 * @see setCoeff(int)
		 * @see clearCoeff(int)
		 */
		inline void setCoeff( int i , bool c ) {
			if ( c ) { setCoeff(i); } else { clearCoeff(i); }
		}

		/**
		 * @brief Adds the specified polynomial to this polynomial and returns
		 *        the result.
		 *
		 * @details The addition of two binary polynomials is essentially the
		 *          same as performing a bit-wise xor operation of its
		 *          representing integers.
		 *
		 * @param f summand
		 *
		 * @return sum of this polynomial and the polynomial <code>f</code>
		 */
		inline SmallBinaryPolynomial add
		( const SmallBinaryPolynomial & f ) const {
			return this->rep^f.rep;
		}

		/**
		 * @brief Subtracts the specified polynomial to this polynomial and
		 *        returns the result.
		 *
		 * @details The subtraction of two binary polynomials is essentially
		 *          the same as performing a bit-wise xor operation of its
		 *          representing integers. In fact, there is no difference
		 *          to addition but for consistency we provide a subtraction
		 *          function.
		 *
		 * @param f subtrahend
		 *
		 * @return difference of this polynomial and the polynomial
		 *         <code>f</code>
		 */
		inline SmallBinaryPolynomial sub
		( const SmallBinaryPolynomial & f ) const {
			return this->rep^f.rep;
		}

		/**
		 * @brief Multiplies this polynomial with the specified polynomial and
		 *        returns the result.
		 *
		 * @details This functions is a wrapper around the carry-less
		 *          multiplication function
		 *          \link MathTools::clmul(uint32_t,uint32_t)\endlink.
		 *          The correctness of the result is thus only guaranteed if
		 *          this polynomial and <code>f</code> are both represented by
		 *          an 32-bit integer.
		 *
		 * @warning If this polynomial or <code>f</code> have degree larger than
		 *          31 then the result of this function may be wrong.
		 *
		 * @param f factor
		 *
		 * @return product of this polynomial and the polynomial
		 *         <code>f</code>
		 */
		inline SmallBinaryPolynomial mul
		( const SmallBinaryPolynomial & f ) const {

			return MathTools::clmul((uint32_t)(this->rep),(uint32_t)(f.rep));
		}

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
		static void divRem
		( SmallBinaryPolynomial & q , SmallBinaryPolynomial & r ,
		  const SmallBinaryPolynomial & a , const SmallBinaryPolynomial & b );

		/**
		 * @brief     Computes the quotient of the division of this
		 *            polynomial divided by <code>f</code>
		 *
		 * @param f
		 *            denominator polynomial
		 *
		 * @return
		 *            quotient of this polynomial divided by <code>f</code>
		 *
		 * @warning
		 *            If <code>f</code> is the zero polynomial an error message
		 *            will be printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 *
		 */
		inline SmallBinaryPolynomial div( const SmallBinaryPolynomial & f ) const {

			SmallBinaryPolynomial q , r;
			divRem(q,r,*this,f);
			return q;
		}

		/**
		 * @brief     Computes the remainder of the division of this
		 *            polynomial divided by <code>f</code>
		 *
		 * @param f
		 *            denominator polynomial
		 *
		 * @return
		 *            remainder of this polynomial modulo <code>f</code>
		 *
		 * @warning
		 *            If <code>f</code> is the zero polynomial an error message
		 *            will be printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 *
		 */
		inline SmallBinaryPolynomial rem
		( const SmallBinaryPolynomial & f ) const {

			SmallBinaryPolynomial q , r;
			divRem(q,r,*this,f);
			return r;
		}

		/**
		 * @brief
		 *            Multiplies this polynomial with <code>f</code> and
		 *            returns the result modulus <code>modulus</code>
		 *
		 * @param f
		 *            Factor polynomial
		 * @param modulus
		 * 		      Modulus polynomial
		 *
		 * @return
		 *            the remainder of this polynomial multiplied with
		 *            <code>f</code> modulo <code>modulus</code>
		 *
		 * @warning
		 *           If <code>modulus</code> is the zero polynomial an
		 *           error message is printed to <code>stderr</code> and the
		 *           program exits with status 'EXIT_FAILURE'.
		 */
		inline SmallBinaryPolynomial mulMod
		( const SmallBinaryPolynomial & f ,
		  const SmallBinaryPolynomial & modulus ) const {

			return mul(f).rem(modulus);
		}

		/**
		 * @brief
		 *          Computes the inverse of this polynomial modulo
		 *          <code>modulus</code> and returns the result.
		 *
		 * @return
		 *          The inverse of this polynomial modulo
		 *          <code>modulus</code>
		 *
		 * @warning
		 *          In the following cases the function prints an error
		 *          message to <code>stderr</code> and exits with status 'EXIT_FAILURE':
		 *          <ul>
		 *           <li>the degree of the modulus is smaller than 1</li>
		 *           <li>this polynomial is the zero polynomial</li>
		 *           <li>this polynomial and the modulus have common divisors
		 *               (except 1)</li>
		 *          </ul>
		 *
		 */
		inline SmallBinaryPolynomial invMod
		( const SmallBinaryPolynomial & modulus ) const {

			if ( modulus.rep < 2 ) {
				std::cerr
				<< "SmallBinaryPolynomial::invMod: Modulus must have degree "
				<< "larger than 0." << std::endl;
				exit(EXIT_FAILURE);
			}

			if ( this->rep == 0 ) {
				std::cerr
				   << "SmallBinaryPolynomial::invMod: Division by zero"
				   << std::endl;
				exit(EXIT_FAILURE);
			}

			SmallBinaryPolynomial g , s , t;

			xgcd(g,s,t,*this,modulus);

			if ( g.rep != 1 ) {
				std::cerr
				<< "SmallBinaryPolynomial::invMod: "
				<< "Polynomial has common divisor with modulus"
				<< std::endl;
				exit(EXIT_FAILURE);
			}

			return s.rem(modulus);
		}

		/**
		 * @brief
		 *            Performs the Extended Euclidean Algorithm.
		 *
		 * @details
		 *            Computes <code>g</code>, <code>s</code>, and
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
		static void xgcd
		( SmallBinaryPolynomial & g ,
		  SmallBinaryPolynomial & s , SmallBinaryPolynomial & t ,
		  const SmallBinaryPolynomial & a , const SmallBinaryPolynomial & b);

		/**
		 * @brief Selects an irreducible polynomial of specified degree
		 *        from a precomputed list.
		 *
		 * @details
		 *        The degree must be between 0 and 31. Otherwise an error
		 *        message is printed to <code>stderr</code> and the program
		 *        exits with status 'EXIT_FAILURE'.
		 *        <br><br>
		 *        The result will be a polynomial that is irreducible and
		 *        that is of the specified degree. The function simply selects
		 *        the polynomial from a fixed list in where for each degree
		 *        (between 1 and 30) there is exactly one polynomial. This
		 *        is in particular useful for computing with finite fields.
		 *
		 * @param degree
		 *             The degree of the sought polynomial
		 *
		 * @return     An irreducible polynomial of specified degree
		 *
		 * @warning
		 *        If <code>degree</code> is not between 1 and 31 an error
		 *        message is printed to <code>stderr</code> and the program
		 *        exits with status 'EXIT_FAILURE'.
		 */
		static SmallBinaryPolynomial irreducible( int degree );

		/**
		 * @brief
		 *            Generates a random polynomial of specified size.
		 *
		 * @details
		 *            The result is a random polynomial of degree smaller
		 *            than <code>size</code>.
		 *
		 * @param size
		 *            Bound of the degree of the random polynomial
		 *
		 * @param tryRandom
		 *            If <code>true</code> there functions attempts to
		 *            use a cryptographically secure random generator (see
		 *            <code>MathTools::rand8()</code>); otherwise
		 *            if <code>false</code> this function wraps around
		 *            the standard random generator <code>rand()</code>s
		 */
		static inline SmallBinaryPolynomial rand
		( int size , bool tryRandom = false ) {

			if ( size < 0 || size > 63) {
				std::cerr
						<< "SmallBinaryPolynomial::rand(int,bool): "
						<< "The size must be between 0 and 63" << std::endl;
				exit(EXIT_FAILURE);
			}

			return MathTools::rand64(tryRandom)%(((uint64_t)1)<<size);
		}
	};

}


#endif /* THIMBLE_SMALLBINARYPOLYNOMIAL_H_ */
