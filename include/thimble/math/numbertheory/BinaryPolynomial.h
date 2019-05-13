/*
 *  THIMBLE --- Research Library for Development and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
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
 * @file BinaryPolynomial.h
 *
 * @brief
 *            Provides a class for representing and computing
 *            with polynomials of arbitrary degree.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_BINARYPOLYNOMIAL_H
#define THIMBLE_BINARYPOLYNOMIAL_H

#include <stdint.h>

#include <fstream>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

    /**
     * @brief
     *           Instances of this class represent binary polynomials
     *           of arbitrary degree.
     */
    class THIMBLE_DLL BinaryPolynomial {

    private:

        /**
         * @brief
         *           Array encoding the coefficients of the polynomial.
         *
         * @details
         *           The array contains \link numWords\endlink unsigned 32 bit
         *           integers data[0],....,data[n-1]. Assume that
         *           data[i]=\f$\sum_{j=0}^{31}b_{i,j}2^j\f$ with
         *           \f$b_{i,j}\in\{0,1\}\f$. Then the polynomial represented
         *           by this instance equals
         *           \f[
         *            \sum_{i=0}^{n-1}X^{i\cdot 32}\sum_{j=0}^{31}b_{i,j}X^j.
         *           \f]
         *           Note that the degree of the polynomial is not necessarily
         *           equals <i>n*32-1</i> as leading coefficients can
         *           be zero.
         */
        uint32_t *data;

        /**
         * @brief
         *           The length of valid unsigned 32 bit integers listed
         *           in the array \link data\endlink.
         *
         * @details
         *           The value of this field should be always at least 1.
         */
        int numWords;

        /**
         * @brief
         *           The degree of the polynomial.
         *
         * @details
         *           Specifies the degree of the polynomial.
         */
        int degree;

    public:

        /**
         * @brief
         *           Standard constructor.
         *
         * @details
         *           If in the degree of the polynomials that this class is used
         *           to represent is known in advance, it can be specified
         *           via the parameter <i>d</i>.
         *
         * @param d
         *           Ensures that the instance can represent polynomials up to
         *           degree <i>d</i> without causing reallocation.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         *
         */
        BinaryPolynomial( int d = 0 );

        /**
         * @brief
         *           Copy constructor.
         *
         * @param f
         *           The binary polynomial of which a copy is constructed.
         *
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        BinaryPolynomial( const BinaryPolynomial & f );

        /**
         * @brief
         *           Destructor.
         *
         * @details
         *           Frees the memory used by this instances to store
         *           the coefficients of the binary polynomial.
         */
        ~BinaryPolynomial();

        /**
         * @brief
         *           Assignment operator (procedural version).
         *
         * @param f
         *           The binary polynomial of which this instances is assigned
         *           a copy.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        void assign( const BinaryPolynomial & f );

        /**
         * @brief
         *           Assignment operator.
         *
         * @param f
         *           The binary polynomial of which this instances is assigned
         *           a copy.
         *
         * @return
         *           A reference to this instance (after assignment).
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryPolynomial & operator=( const BinaryPolynomial & f ) {
            assign(f);
            return *this;
        }

        /**
         * @brief
         *           Accesses the degree of the binary polynomial.
         *
         * @return
         *           The degree of the binary polynomial.
         */
        inline int deg() const {
            return this->degree;
        }

        /**
         * @brief
         *           Returns the Hamming weight of this polynomial.
         *
         * @details
         *           Write the polynomial as
         *           \f[
         *            \sum_{j=0}^d c_j\cdot X^j
         *           \f]
         *           where \f$c_j\in\{0,1\}\f$. The function returns
         *           the number of \f$j\f$ where \f$c_j\neq 0\f$. More
         *           specifically, it returns the cardinality
         *           \f[
         *            \#\{j=0,...,d~|~c_j\neq 0\}.
         *           \f]
         *
         * @return
         *           The Hamming weight of this polynomial.
         */
        int hammingWeight() const;

        /**
         * @brief
         *           Tests whether the polynomial is identically 0.
         *
         * @return
         *           <code>true</code> if the polynomial is zero and
         *           <code>false</code> otherwise.
         */
        inline bool isZero() const {
            return this->degree < 0;
        }

        /**
         * @brief
         *           Clears all data held by this polynomial
         *           such that it becomes identically to 0.
         */
        void setZero();

        /**
         * @brief
         *           Tests whether the polynomial is identically 1.
         *
         * @return
         *           <code>true</code> if the polynomial equals 1 and
         *           <code>false</code> otherwise.
         */
        inline bool isOne() const {
            return this->degree == 0;
        }

        /**
         * @brief
         *           Assignment of this instance to be the constant
         *           polynomial equals 1.
         */
        inline void setOne() {
            setZero();
            this->data[0] = 1;
            this->degree = 0;
        }

        /**
         * @brief
         *           Tests whether this polynomial represents the
         *           monomial <i>X</i>.
         *
         * @details
         *           <code>true</code> if this polynomial equals the
         *           monomial <i>X</i> and otherwise <code>false</code>.
         */
        inline bool isX() const {
            return this->degree == 1 && this->data[0] == 2;
        }

        /**
         * @brief
         *           Assigns the monomial <i>X</i> to this polynomial.
         */
        inline void setX() {
            setZero();
            this->data[0] = 2;
            this->degree = 1;
        }

        /**
         * @brief
         *           Access the <i>i</i>th coefficients of the binary
         *           polynomial.
         *
         * @details
         *           Write the polynomial as
         *           \f[
         *            \sum_{j=0}^dc_j\cdot X^j
         *           \f]
         *           with \f$c_j\in\{0,1\}\f$. The function returns
         *           the coefficient \f$c_i\f$.
         *
         * @param i
         *           The index of the coefficient that is accessed.
         *
         * @return
         *           <code>true</code> if \f$c_i=1\f$; otherwise, if
         *           \f$c_i=0\f$, the function returns <code>false</code>.
         *
         * @warning
         *           If <i>i</i> is less than 0, an error message is
         *           printed to <code>stderr</code> and the program
         *           exits with status 'EXIT_FAILURE'.
         */
        bool getCoeff( int i ) const;

        /**
         * @brief
         *           Ensures that the <i>i</i>th coefficient of the
         *           polynomial is set to 1.
         *
         * @details
         *           Write the polynomial as
         *           \f[
         *            \sum_{j=0}^dc_j\cdot X^j
         *           \f]
         *           with \f$c_j\in\{0,1\}\f$. The method ensures
         *           that \f$c_i=1\f$.
         *
         * @param i
         *           The index of the coefficient to ensuring to be
         *           equal 1.
         *
         * @warning
         *           If <i>i</i> is less than 0, an error message is
         *           printed to <code>stderr</code> and the program
         *           exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        void setCoeff( int i );

        /**
         * @brief
         *           Ensures that the <i>i</i>th coefficient of the
         *           polynomial is set to 0.
         *
         * @details
         *           Write the polynomial as
         *           \f[
         *            \sum_{j=0}^dc_j\cdot X^j
         *           \f]
         *           with \f$c_j\in\{0,1\}\f$. The method ensures
         *           that \f$c_i=0\f$.
         *
         * @param i
         *           The index of the coefficient to ensuring to be
         *           equal 1.
         *
         * @warning
         *           If <i>i</i> is less than 0, an error message is
         *           printed to <code>stderr</code> and the program
         *           exits with status 'EXIT_FAILURE'.
         */
        void clearCoeff( int i );

        /**
         * @brief
         *           Updates the <i>i</i>th coefficient by the
         *           specified value.
         * @details
         *           Write the polynomial as
         *           \f[
         *            \sum_{j=0}^dc_j\cdot X^j
         *           \f]
         *           with \f$c_j\in\{0,1\}\f$. The method updates
         *           the polynomial such that it represents
         *           \f[
         *            c\cdot X^i+\sum_{j=0,j\neq i}^dc_j\cdot X^j.
         *           \f]
         *
         * @param i
         *           The index of the coefficient of the polynomial
         *           being updated.
         *
         * @param c
         *           Specifies the value of the <i>i</i>th coefficient
         *           of the updated polynomial which is <code>TRUE</code>
         *           or <code>FALSE</code> if the specified coefficient
         *           is 1 or 0, respectively.
         *
         * @warning
         *           If <i>i</i> is less than 0, an error message is
         *           printed to <code>stderr</code> and the program
         *           exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline void setCoeff( int i , bool c ) {

            if ( c ) {
                setCoeff(i);
            } else {
                clearCoeff(i);
            }
        }

        /**
         * @brief
         *           Ensures that this polynomial can store at least
         *           <i>d+1</i> coefficients without reallocation.
         *
         * @param d
         *           A degree bound for a polynomial that can
         *           be represented by this instance without
         *           reallocation.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        void ensureDegree( int d );

        /**
         * @brief
         *           Initializes this instance by a polynmial with the
         *           specified number of random bit coefficients.
         *
         * @details
         *           The degree of the resulting polynomial will be
         *           strictly smaller than <code>numCoefficients</code>.
         *
         * @param numCoefficients
         *           The number of random coefficients.
         *
         * @param tryRandom
         *           If <code>true</code>, the method uses a cryptographic
         *           number generator if available on the system; otherwise,
         *           the method wraps around the standard <code>rand()</code>
         *           function.
         *
         * @warning
         *           If <code>numCoefficients</code> is negative, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        void random( int numCoefficients  , bool tryRandom = false );

        /**
         * @brief
         *           Replaces this instance by polynomial with the
         *           specified number of random coefficients such that
         *           it has the specified Hamming weight.
         *
         * @param numCoefficients
         *           The number of random coefficients.
         *
         * @param hammingWeight
         *           The Hamming weight of the random polynomial.
         *
         * @param tryRandom
         *           If <code>true</code>, the method uses a cryptographic
         *           number generator if available on the system; otherwise,
         *           the method wraps around the standard <code>rand()</code>
         *           function.
         *
         * @warning
         *           If <code>hammingWeight</code> is negative or greater than
         *           <code>numCoefficients</code>, an error message is printed
         *           to <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If <code>numCoefficients</code> is negative, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        void wrandom
        ( int numCoefficients , int hammingWeight ,
          bool tryRandom = false );

        /**
         * @brief
         *           Replaces this polynomial by its multipliciation
         *           with \f$X^e\f$.
         *
         * @details
         *           Write the polynomial as
         *           \f[
         *            \sum_{j=0}^dc_j\cdot X^j
         *           \f]
         *           with \f$c_j\in\{0,1\}\f$. The method replaces the
         *           polynomial by
         *           \f[
         *            X^e\cdot \sum_{j=0}^dc_j\cdot X^j.
         *           \f]
         *
         * @param e
         *           The exponent of the shift.
         *
         * @warning
         *           If <i>e</i> is negative, an error message is printed
         *           to <code>stderr</code> and the program exits with
         *           status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        void leftShift( int e );

        /**
         * @brief
         *           Essentially, polynomial division by the
         *           monomial \f$X^e\f$.
         *
         * @details
         *           Write the polynomial as
         *           \f[
         *            \sum_{j=0}^dc_j\cdot X^j
         *           \f]
         *           with \f$c_j\in\{0,1\}\f$. The method replaces the
         *           polynomial by
         *           \f[
         *            \sum_{j=e}^dc_j\cdot X^{j-e}.
         *           \f]
         *
         * @param e
         *           The exponent of the shift.
         *
         * @warning
         *           If <i>e</i> is negative, an error message is printed
         *           to <code>stderr</code> and the program exits with
         *           status 'EXIT_FAILURE'.
         */
        void rightShift( int e );

        /**
         * @brief
         *           Replaces the polynomial by its derivative.
         *
         * @details
         *           Write the polynomial as
         *           \f[
         *            \sum_{j=0}^d c_j\cdot X^j
         *           \f]
         *           with \f$c_j\in\{0,1\}\f$. The method replaces
         *           the polynomial by the formal derivative
         *           \f[
         *            \sum_{j=1,3,5,...,d} c_j\cdot X^{j-1}.
         *           \f]
         */
        void derive();

        /**
         * @brief
         *           Computes the <i>d</i>-reverse of this polynomial.
         *
         * @details
         *           The <i>d</i>-reversal of a polynomial
         *           \f[
         *            \sum_{j=0}^n f_j\cdot X^j
         *           \f]
         *           is defined as the polynomial
         *           \f[
         *            \sum_{j=0}^n f_j\cdot X^{d-j}.
         *           \f]
         *
         * @param d
         *           Positive integer greater than or equals the polynomial's
         *           degree.
         *
         * @return
         *           The <i>d</i>-reverse of this polynomial.
         *
         * @warning
         *           If <i>d</i> is negative or smaller than the
         *           polynomial's degree (i.e., \link deg()\endlink), an
         *           error message is printed to <code>stderr</code> and
         *           the program exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        BinaryPolynomial reverse( int d ) const;

        /**
         * @brief
         *            Determines the <i>n</i>th cyclotomic polynomial.
         *
         * @details
         *            The <i>n</i>th cyclotomic polynomial is defined as
         *            the following:
         *            \f[
         *   \Phi_n(X):=\prod_{1\leq k<n,\\\gcd(k,n)=1}(X-e^{2\pi ik/n}).
         *            \f]
         *            This function computes this polynomial whose
         *            coefficients are reduced modulo 2.
         *
         *            The implementation of the cyclotomic polynomial
         *            computation follows the description of
         *            <ul>
         *             <li><b>[vzGth]</b> v.z. Gathen and Gerhard (2003).
         *              <i>Modern Computer Algebra</i>. Cambridge University
         *              Press, Cambridge (UK), 2nd edition.
         *             </li>
         *            </ul>
         *
         * @param n
         *            Index of the cyclotomic polynomial.
         *
         * @return
         *            The <i>n</i>th cyclotomic polynomial.
         *
         * @warning
         *            If <i>n</i> is zero or negative, an error message is
         *            printed to <code>stderr</code> and the program exits
         *            with status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static BinaryPolynomial cyclo( int n );


        /**
         * @brief
         *           Tests two polynomials on equality.
         *
         * @param f
         *           First input polynomial.
         *
         * @param g
         *           Second input polynomial.
         *
         * @return
         *           <code>true</code> if <i>f</i> and <i>g</i> represent
         *           the same polynomials and, otherwise, <code>false</code>.
         */
        static bool areEqual
        ( const BinaryPolynomial & f , const BinaryPolynomial & g );

        /**
         * @brief
         *           Determines the number of coefficients in which two
         *           polynomials differ.
         *
         * @details
         *           Write
         *           \f[
         *            f(X)=\sum_{j=0}^n f_j X^j
         *           \f]
         *           and
         *           \f[
         *            g(X)=\sum_{j=0}^n g_j X^j
         *           \f]
         *           where \f$f_j,g_j\in\{0,1\}\f$. The function returns
         *           the number of <i>j</i> in which the \f$f_j\f$ and
         *           \f$g_j\f$ differ; more specifically, it returns
         *           \f[
         *            \#\{j=0,...,n~|~f_j\neq g_j\}.
         *           \f]
         *
         * @param f
         *           First input polynomial.
         *
         * @param g
         *           Second input polynomial.
         *
         * @return
         *           The hamming distance between <i>f</i> and <i>g</i>.
         */
        static int hammingDistance
        ( const BinaryPolynomial & f , const BinaryPolynomial & g );

        /**
         * @brief
         *           Swaps the content of two polynomials such the
         *           the one represents the other.
         *
         * @param f
         *           Polynomial being assigned the fields of <i>g</i>.
         *
         * @param g
         *           Polynomial being assigned the fields of <i>f</i>.
         */
        static void swap( BinaryPolynomial & f , BinaryPolynomial & g );

        /**
         * @brief
         *           Computes the sum of two binary polynomials.
         *
         * @details
         *           The method replaces the polynomial <i>f</i>
         *           by the sum of the two polynomials <i>g</i>
         *           and <i>h</i>.
         *
         * @param f
         *           Sum of <i>g</i> and <i>h</i>.
         *
         * @param g
         *           First summand.
         *
         * @param h
         *           Second summand.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void add
        ( BinaryPolynomial & f ,
          const BinaryPolynomial & g ,
          const BinaryPolynomial & h );

        /**
         * @brief
         *           Computes the difference of two binary polynomials.
         *
         * @details
         *           The method replaces the polynomial <i>f</i>
         *           by the difference of the two polynomials <i>g</i>
         *           and <i>h</i>.
         *
         *           Note that, since we are computing with binary
         *           arithmetic, computating the difference between
         *           two polynomials is equivalent to addition.
         *
         * @param f
         *           Difference of <i>g</i> and <i>h</i>.
         *
         * @param g
         *           Minuend.
         *
         * @param h
         *           Subtrahend.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline static void sub
        ( BinaryPolynomial & f ,
          const BinaryPolynomial & g ,
          const BinaryPolynomial & h ) {
            add(f,g,h);
        }

        /**
         * @brief
         *           Computes the produce of two binary polynomials.
         *
         * @details
         *           The method replaces the polynomial <i>f</i>
         *           by the product of the two polynomials <i>g</i>
         *           and <i>h</i>.
         *
         * @param f
         *           Product of <i>g</i> and <i>h</i>.
         *
         * @param g
         *           First factor.
         *
         * @param h
         *           Second factor.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void mul
        ( BinaryPolynomial & f ,
          const BinaryPolynomial & g ,
          const BinaryPolynomial & h );

        /**
         * @brief
         *           Computes a Euclidean division with remainder.
         *
         * @details
         *           Given two polynomials <i>a</i> and a non-zero <i>b</i>
         *           there exits unique polynomials <i>q</i> and <i>r</i>
         *           such that
         *           \f$
         *            a=q\cdot b+r
         *           \f$
         *           with \f$\deg(r)<\deg(b)\f$. These two polynomials
         *           are determined by the method.
         *
         * @param q
         *           The quotient of the Euclidean division.
         *
         * @param r
         *           The remainder of the Euclidean division.
         *
         * @param a
         *           Numerator.
         *
         * @param b
         *           Denominator.
         *
         * @warning
         *           If <i>b</i> is zero, an error message is printed
         *           to <code>stderr</code> and the program exits with
         *           status 'EXIT_FAILURE'.
         *
         * @warning
         *           If <i>q</i> and <i>r</i> are the same references,
         *           an error message is printed to <code>stderr</code>
         *           and the program exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void divRem
        ( BinaryPolynomial & q , BinaryPolynomial & r ,
          const BinaryPolynomial & a , const BinaryPolynomial & b );

        /**
         * @brief
         *           Computes the inverse of a polynomial modulo
         *           another polynomial.
         *
         * @details
         *           The method attempts to find a polynomial <i>s</i>
         *           of degree smaller than \f$\deg(b)\f$ such that
         *           \f[
         *            1=c(X)\cdot a(X)~rem~b.
         *           \f]
         *           If such a polynomial exists, the function
         *           returns <code>true</code> and replaces the content
         *           of <i>c</i> correspondingly; otherwise, if the inverse
         *           does not exist, the function returns <code>false</code>
         *           and leaves the content of <i>c</i> unchanged.
         *
         * @param c
         *           Output of the inverse of <i>a</i> modulo <i>b</i>.
         *
         * @param a
         *           Input polynomial of which the inverse is computed.
         *
         * @param b
         *           Input modulus polynomial.
         *
         * @return
         *           <code>true</code> if the inverse have been determined
         *           successfully, in particular, if it exists; otherwise
         *           the function returns <code>false</code>.
         *
         * @warning
         *           If <i>a</i> and <i>b</i> are both zero, the program
         *           prints an error message to <code>stderr</code> and
         *           exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static bool invMod
        ( BinaryPolynomial & c ,
          const BinaryPolynomial & a , const BinaryPolynomial & b );

        /**
         * @brief
         *           Computes the quotient of two polynomials via Euclidean
         *           division.
         *
         * @details
         *           The method calls the \link divRem()\endlink
         *           method and dismisses the remainder.
         *
         * @param q
         *           Quotient.
         *
         * @param a
         *           Numerator.
         *
         * @param b
         *           Denominator.
         *
         * @see divRem()
         *
         * @warning
         *           If <i>b</i> is zero, an error message is printed
         *           to <code>stderr</code> and the program exits with
         *           status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline static void div
        ( BinaryPolynomial & q , const BinaryPolynomial & a ,
          const BinaryPolynomial & b ) {

            BinaryPolynomial tmp;
            divRem(q,tmp,a,b);
        }

        /**
         * @brief
         *           Computes the remainder of two polynomials via Euclidean
         *           division.
         *
         * @details
         *           The method calls the \link divRem()\endlink
         *           method and dismisses the quotient.
         *
         * @param r
         *           Remainder.
         *
         * @param a
         *           Numerator.
         *
         * @param b
         *           Denominator.
         *
         * @see divRem()
         *
         * @warning
         *           If <i>b</i> is zero, an error message is printed
         *           to <code>stderr</code> and the program exits with
         *           status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline static void rem
        ( BinaryPolynomial & r , const BinaryPolynomial & a ,
          const BinaryPolynomial & b ) {

            BinaryPolynomial tmp;
            divRem(tmp,r,a,b);
        }

        /**
         * @brief
         *          Computes the greatest common divisor of two binary
         *          polynomials.
         *
         * @param g
         *          The greatest common divisors of <i>a</i> and <i>b</i>.
         *
         * @param a
         *          First polynomial.
         *
         * @param b
         *          Second polynomial.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void gcd
        ( BinaryPolynomial & g ,
          const BinaryPolynomial & a , const BinaryPolynomial & b );

        /**
         * @brief
         *           Performs the extended Euclidean algorithm.
         *
         * @details
         *           The method determines the greatest common divisor
         *           <i>g</i> of the two polynomials <i>a</i> and <i>b</i>
         *           and co-prime polynomials <i>s</i> and <i>t</i> such that
         *           \f[
         *            g=s\cdot a+t\cdot b.
         *           \f]
         *
         * @param g
         *           Output of the greatest common divisor of <i>a</i>
         *           and <i>b</i>.
         *
         * @param s
         *           Output polynomial as described above.
         *
         * @param t
         *           Output polynomial as described above.
         *
         * @param a
         *           First input polynomial.
         *
         * @param b
         *           Second input polynomial.
         *
         * @warning
         *           If <i>g</i>, <i>s</i> or <i>t</i> are of the same
         *           reference, an error message is printed to
         *           <code>stderr</code> and the program exits with
         *           status 'EXIT_FAILURE'.
         *
         * @warning
         *           If <i>a</i> and <i>b</i> are both zero, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void xgcd
        ( BinaryPolynomial & g ,
          BinaryPolynomial & s , BinaryPolynomial & t ,
          const BinaryPolynomial & a , const BinaryPolynomial & b );

        /**
         * @brief
         *            Computes the partial greatest common divisor of two
         *            polynomials.
         *
         * @details
         *            This methods performs the traditional extended euclidean
         *            algorithm, i.e., determines \f$s,t\f$ such that
         *            \f$g=s\cdot a+t\cdot b\f$, until the degree of
         *            \f$g\f$ becomes smaller than \f$d\f$.
         *
         * @param g
         *            Will contain the partial greatest common divisor
         *            of <i>a</i> and <i>b</i>.
         *
         * @param s
         *            see details
         * @param t
         *            see details
         * @param a
         *            see details
         * @param b
         *            see details
         * @param d
         *            see details
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void pgcd
        ( BinaryPolynomial & g ,
          BinaryPolynomial & s ,
          BinaryPolynomial & t ,
          const BinaryPolynomial & a ,
          const BinaryPolynomial & b ,
          int d );

        /**
         * @brief
         *            Computes the <i>(n,n-k)</i>-Pade approximant to the
         *            polynomial.
         *
         * @details
         *            The method determines co-prime polynomials <i>r</i> and
         *            <i>t</i> with \f$\deg(r)<k\f$ and \f$\deg(t)\leq n-k\f$
         *            such that
         *            \f[
         *             \frac{r(X)}{t(X)}=g(X)~mod~X^n
         *            \f]
         *            where \f$t(X)\f$ is not divisible by \f$X\f$.
         *
         *            For more details on <em>Pade approximation</em> we refer
         *            to Section 5.9 in
         *            <ul>
         *             <li><b>[vzGth]</b> v.z. Gathen and Gerhard (2003).
         *              <i>Modern Computer Algebra</i>. Cambridge University
         *              Press, Cambridge (UK), 2nd edition.
         *             </li>
         *            </ul>
         *
         * @param r
         *            Output numerator polynomial.
         *
         * @param t
         *            Output denominator polynomial.
         *
         * @param g
         *            Input polynomial.
         *
         * @param n
         *            Bounds the degree of <i>r</i>.
         *
         * @param k
         *            Bounds the degree of <i>t</i>, i.e., such that
         *            \f$\deg(t)\leq n-k\f$.
         *
         * @warning
         *            If <i>n</i> or <i>k</i> is negative, the degree
         *            of <i>g</i> is greater than or equals <i>n</i>, or
         *            if <i>k>n</i>, an error message is written to
         *            <code>stderr</code> and the program exits with
         *            status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void pade
        ( BinaryPolynomial & r , BinaryPolynomial & t ,
          const BinaryPolynomial & g , int n , int k );

        /**
         * @brief
         *           Computes the minimal polynomial of a polynomial modulo
         *           another polynomial.
         *
         * @details
         *           The minimal polynomial of <i>g</i> modulo <i>f</i> is
         *           defined as the non-zero polynomial <i>h</i> of minimal
         *           degree such that
         *           \f[
         *            h(g(X))~mod~f(X)=0.
         *           \f]
         *
         * @param h
         *           Output to contain the minimal polynomial of <i>g</i>
         *           modulo <i>f</i>.
         *
         * @param g
         *           Input polynomial.
         *
         * @param f
         *           Input modulus polynomial; must be of degree greater
         *           than 1.
         *
         * @warning
         *           If <i>f</i> is of degree smaller than 1, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void minPolyMod
        ( BinaryPolynomial & h ,
          const BinaryPolynomial & g , const BinaryPolynomial & f );

        /**
         * @brief
         *          Computes the composition of two binary polynomials.
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
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void eval
        ( BinaryPolynomial & h ,
          const BinaryPolynomial & f ,
          const BinaryPolynomial & g );

        /**
         * @brief
         *           Computes the composition of two polynomials being
         *           reduced by a modulus.
         *
         * @details
         *           The method performs <i>h(X)=f(g(X))</i> mod <i>m</i>.
         *
         * @param h
         *           Modular composition.
         *
         * @param f
         *           The outer polynomial.
         *
         * @param g
         *           The inner polynomial.
         *
         * @param m
         *           Modulus.
         *
         * @warning
         *           If <i>m</i> is zero, an error message is printed
         *           to <code>stderr</code> and the program exits with
         *           status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void evalMod
        ( BinaryPolynomial & h , const BinaryPolynomial & f ,
          const BinaryPolynomial & g , const BinaryPolynomial & m );

        /**
         * @brief
         *          Computes the square root of a polynomial.
         *
         * @details
         *          If the input polynomial <i>f</i> is the square
         *          of another polynomial <i>g</i>, i.e., \f$f=g^2\f$, then
         *          the method will correctly determined <i>g</i>.
         *
         * @param g
         *          On output square root of the input polynomial <i>f</i>.
         *
         * @param f
         *          On input the square of the polynomial <i>g</i>.
         *
         * @see isSquare()
         *
         * @warning
         *          If <i>f</i> is not the square of a polynomial, the
         *          result of the method is undocumented.
         *
         * @warning
         *          If not enough memory could be provided, an error
         *          message is printed to <code>stderr</code> and the
         *          program exits with status 'EXIT_FAILURE'.
         */
        static void squareRoot
            ( BinaryPolynomial & g , const BinaryPolynomial & f );

        /**
         * @brief
         *          Tests whether this polynomial is the square of
         *          another binary polynomial.
         *
         * @return
         *          <code>true</code> if this polynomial is a square
         *          of another polynomial and <code>false</code> otherwise.
         */
        bool isSquare() const;

        /**
         * @brief
         *          Determines whether the binary polynomial is irreducible.
         *
         * @details
         *          The implementation of the irreducible test follows
         *          the description of Algorithm 14.36 in
         *            <ul>
         *             <li><b>[vzGth]</b> v.z. Gathen and Gerhard (2003).
         *              <i>Modern Computer Algebra</i>. Cambridge University
         *              Press, Cambridge (UK), 2nd edition.
         *             </li>
         *            </ul>
         *
         * @return
         *          <code>true</code> if the polynomial is irreducible and
         *          otherwise, if not irreducible, <code>false</code>.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        bool isIrreducible() const;

        /**
         * @brief
         *           Performs <i>Equal Degree Splitting</i> to split
         *           a proper factor from the polynomial as it is
         *           known to be a product of irreducible distinct
         *           polynomials of degree exactly <i>d</i>.
         *
         * @details
         *           For details of the implemented algorithm we refer
         *           to Algorithm 14.8 and Excercise 14.16 in
         *            <ul>
         *             <li><b>[vzGth]</b> v.z. Gathen and Gerhard (2003).
         *              <i>Modern Computer Algebra</i>. Cambridge University
         *              Press, Cambridge (UK), 2nd edition.
         *             </li>
         *            </ul>
         *
         * @param d
         *           The degree of all irreducible divisors of the polynomial.
         *
         * @return
         *           A proper divisor of the polynomial (of which degree
         *           is a multiple of <i>d</i>).
         *
         * @warning
         *           If the polynomial is not the product of distinct
         *           irreducible polynomials of degree <i>d</i>, the
         *           function runs into undocumented behaviour.
         *
         * @warning
         *           If <i>d</i> is zero, negative or equals the degree
         *           of this polynomial or if <i>d</i> does not divide
         *           the degree, an error message is printed
         *           <code>stderr</code> an the program exits with
         *           status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        BinaryPolynomial equalDegreeSplitting( int d ) const;

        /**
         * @brief
         *           Determines an irreducible divisors from the polynomial
         *           for which it is know that it is a product of irreducible
         *           and distinct polynomials of degree exactly <i>d</i>.
         *
         * @param d
         *           The degree of all irreducible divisors of the polynomial.
         *
         * @return
         *           A proper divisor of the polynomial of degree
         *           <i>d</i> (which is assumed to be irreducible).
         *
         * @warning
         *           If the polynomial is not the product of distinct
         *           irreducible polynomials of degree <i>d</i>, the
         *           function runs into undocumented behaviour.
         *
         * @warning
         *           If <i>d</i> is zero, negative or does not divide
         *           the degree of this polynomial, an error message is
         *           printed <code>stderr</code> an the program exits with
         *           status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        BinaryPolynomial splitIrreducibleEDF( int d ) const;

        /**
         * @brief
         *           Determines an irreducible divisor from the
         *           polynomial which is assumed to be square-free.
         *
         * @return
         *           An irreducible divisor of this polynomial.
         *
         * @warning
         *           If this polynomial is not square-free, the function
         *           runs into undocumented behaviour.
         *
         * @warning
         *           If the degree of this polynomial is smaller than
         *           1, an error message is printed to <code>stderr</code>
         *           an the program exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        BinaryPolynomial splitIrreducibleDDF() const;

        /**
         * @brief
         *           Determines the product of all distinct irreducible
         *           divisors of this polynomial.
         *
         * @return
         *           The product of all distinct and irreducible divisors
         *           of this polynomial.
         *
         * @warning
         *           If the polynomial is zero, an error message is printed
         *           to <code>stderr</code> and the program exits with
         *           status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        BinaryPolynomial squareFreeDecomposition() const;

        /**
         * @brief
         *           Determines an irreducible factor of this polynomial.
         *
         * @details
         *           If the polynomial is irreducible, the returned polynomial
         *           will be equals to this polynomial. This method is useful
         *           if a single irreducible divisor of a polynomial is needed.
         *           For irreducibility tests it is more efficient to use
         *           the \link isIrreducible()\endlink function.
         *
         * @return
         *           An irreducible factor of this polynomial.
         *
         * @see isIrreducible()
         *
         * @warning
         *           If the polynomial is zero, an error message is printed
         *           to <code>stderr</code> and the program exits with
         *           status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        BinaryPolynomial splitIrreducible() const;

        /**
         * @brief
         *           Enables to compare two binary polynomials via the
         *           '=='-operator.
         *
         * @param g
         *           Binary polynomial with which this polynomial is
         *           compared.
         *
         * @return
         *           <code>true</code> if this polynomial equals <i>g</i>
         *           or, otherwise, <code>false</code>.
         */
        inline bool operator==( const BinaryPolynomial & g ) const {
            return areEqual(*this,g);
        }

        /**
         * @brief
         *           Enables to compare two binary polynomials via the
         *           '!='-operator.
         *
         * @param g
         *           Binary polynomial with which this polynomial is
         *           compared.
         *
         * @return
         *           <code>false</code> if this polynomial equals <i>g</i>
         *           or, otherwise, <code>true</code>.
         */
        inline bool operator!=( const BinaryPolynomial & g ) const {
            return !areEqual(*this,g);
        }

        /**
         * @brief
         *           Enables to add two binary polynomials via the
         *           '+'-operator.
         *
         * @param g
         *           Summand.
         *
         * @return
         *           The sum of this polynomial and <i>g</i>.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryPolynomial operator+( const BinaryPolynomial & g ) const {
            BinaryPolynomial h;
            add(h,*this,g);
            return h;
        }

        /**
         * @brief
         *           Enables to compute the difference between two binary
         *           polynomials via the '-'-operator.
         *
         * @param g
         *           Subtrahend.
         *
         * @return
         *           The difference between this polynomial and <i>g</i>.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryPolynomial operator-( const BinaryPolynomial & g ) const {
            BinaryPolynomial h;
            sub(h,*this,g);
            return h;
        }

        /**
         * @brief
         *           Enables to compute the product of two binary polynomials
         *           via the '*'-operator.
         *
         * @param g
         *           Factor.
         *
         * @return
         *           The product of this polynomial and <i>g</i>.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryPolynomial operator*( const BinaryPolynomial & g ) const {
            BinaryPolynomial h;
            mul(h,*this,g);
            return h;
        }

        /**
         * @brief
         *           Enables to compute the quotient of a binary polynomial
         *           divided by another polynomial via the '/'-operator.
         *
         * @param g
         *           Denominator.
         *
         * @return
         *           The quotient of this polynomial divided by <i>g</i>.
         *
         * @warning
         *           If <i>g</i> is zero, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryPolynomial operator/( const BinaryPolynomial & g ) const {
            BinaryPolynomial h;
            div(h,*this,g);
            return h;
        }

        /**
         * @brief
         *           Enables to compute the remainder of a binary polynomial
         *           when divided by another polynomial via the '%'-operator.
         *
         * @param g
         *           Denominator.
         *
         * @return
         *           The remainder of this polynomial divided by <i>g</i>.
         *
         * @warning
         *           If <i>g</i> is zero, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryPolynomial operator%( const BinaryPolynomial & g ) const {
            BinaryPolynomial h;
            rem(h,*this,g);
            return h;
        }

        /**
         * @brief
         *           Enables that a polynomial can be added to another
         *           polynomial via the '+=' operator.
         *
         * @param g
         *           Summand.
         *
         * @return
         *           A reference to this polynomial to which the polynomial
         *           <i>g</i> has been added.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryPolynomial & operator+=( const BinaryPolynomial & g ) {
            BinaryPolynomial::add(*this,*this,g);
            return *this;
        }

        /**
         * @brief
         *           Enables that a polynomial can be subtracted from another
         *           polynomial via the '-=' operator.
         *
         * @param g
         *           Subtrahend.
         *
         * @return
         *           A reference to this polynomial from which the polynomial
         *           <i>g</i> has been subtracted.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryPolynomial & operator-=( const BinaryPolynomial & g ) {
            BinaryPolynomial::sub(*this,*this,g);
            return *this;
        }

        /**
         * @brief
         *           Enables that a polynomial can be multiplied to another
         *           polynomial via the '*=' operator.
         *
         * @param g
         *           Factor.
         *
         * @return
         *           A reference to this polynomial to which the polynomial
         *           <i>g</i> has been multiplied.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryPolynomial & operator*=( const BinaryPolynomial & g ) {
            BinaryPolynomial::mul(*this,*this,g);
            return *this;
        }

        /**
         * @brief
         *           Enables that a polynomial can be divided by another
         *           polynomial via the '/=' operator.
         *
         * @param g
         *           Denominator.
         *
         * @return
         *           A reference to this polynomial which has been replaced
         *           by the quotient of the polynomial represented by this
         *           instance on input divided by <i>g</i>.
         *
         * @warning
         *           If <i>g</i> is zero, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryPolynomial & operator/=( const BinaryPolynomial & g ) {
            BinaryPolynomial::div(*this,*this,g);
            return *this;
        }

        /**
         * @brief
         *           Enables that the remainder of a polynomial modulo
         *           another polynomial can be computed via the '%='-operator.
         *
         * @param g
         *           Denominator.
         *
         * @return
         *           A reference to this polynomial which has been replaced
         *           by the remainder of the polynomial represented by this
         *           instance on input modulo <i>g</i>.
         *
         * @warning
         *           If <i>g</i> is zero, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryPolynomial & operator%=( const BinaryPolynomial & g ) {
            BinaryPolynomial::rem(*this,*this,g);
            return *this;
        }

        /**
         * @brief
         *           Returns the array encoding the coefficients of this
         *           polynomial.
         *
         * @details
         *           The array contains \link deg()\endlink+1 unsigned 32 bit
         *           integers data[0],....,data[n-1]. Assume that
         *           data[i]=\f$\sum_{j=0}^{31}b_{i,j}2^j\f$ with
         *           \f$b_{i,j}\in\{0,1\}\f$. Then the polynomial represented
         *           by this instance equals
         *           \f[
         *            \sum_{i=0}^{n-1}X^{i\cdot 32}\sum_{j=0}^{31}b_{i,j}X^j.
         *           \f]
         *
         * @return
         *           Array encoding the coefficients of this polynomial.
         */
        inline const uint32_t *getData() const {

            return this->data;
        }

        /**
         * @brief
         *           Returns the array encoding the coefficients of this
         *           polynomial.
         *
         * @details
         *            The result of this function is the same as the result of
         *            \link getData()\endlink except that the content of the
         *            polynomial can be changed by changing the content of the
         *            result.
         *
         * @warning
         *            Use the result with caution and perform only changes if
         *            you really know what your are doing; otherwise, you might
         *            experience undocumented behaviour.
         *
         * @return
         *           Array encoding the coefficients of this polynomial.
         */
        inline uint32_t *getData_nonconst() {

            return this->data;
        }

        /**
         * @brief
         *           Updates the degree of the polynomial.
         *
         * @details
         *           After the coefficients stored in
         *           \link getData()_nonconst\endlink
         *           change, the \link deg() degree\endlink of the polynomial
         *           can smaller than indicated. The method updates the degree
         *           correspondingly.
         *
         * @attention
         *           You probably may not want to use this method unless you
         *           have changed coefficients of this polynomial in
         *           combination with the \link getData_nonconst()\endlink
         *           function. The \link getData_nonconst()\endlink function,
         *           however, should only be used if you are really knowing
         *           what you are doing unless you want experience
         *           undocumented behaviour.
         */
        void degreeCouldBeSmaller();
    };

    /**
     * @brief
     *           Tests two polynomials on equality.
     *
     * @param f
     *           First input polynomial.
     *
     * @param g
     *           Second input polynomial.
     *
     * @return
     *           <code>true</code> if <i>f</i> and <i>g</i> represent
     *           the same polynomials and, otherwise, <code>false</code>.
     */
    inline bool areEqual
    ( const BinaryPolynomial & f , const BinaryPolynomial & g ) {

        return BinaryPolynomial::areEqual(f,g);
    }

    /**
     * @brief
     *           Determines the number of coefficients in which two
     *           polynomials differ.
     *
     * @details
     *           Write
     *           \f[
     *            f(X)=\sum_{j=0}^n f_j X^j
     *           \f]
     *           and
     *           \f[
     *            g(X)=\sum_{j=0}^n g_j X^j
     *           \f]
     *           where \f$f_j,g_j\in\{0,1\}\f$. The function returns
     *           the number of <i>j</i> in which the \f$f_j\f$ and
     *           \f$g_j\f$ differ; more specifically, it returns
     *           \f[
     *            \#\{j=0,...,n~|~f_j\neq g_j\}.
     *           \f]
     *
     * @param f
     *           First input polynomial.
     *
     * @param g
     *           Second input polynomial.
     *
     * @return
     *           The hamming distance between <i>f</i> and <i>g</i>.
     */
    inline int hammingDistance
    ( const BinaryPolynomial & f , const BinaryPolynomial & g ) {

        return BinaryPolynomial::hammingDistance(f,g);
    }

    /**
     * @brief
     *           Swaps the content of two polynomials such the
     *           the one represents the other.
     *
     * @param f
     *           Polynomial being assigned the fields of <i>g</i>.
     *
     * @param g
     *           Polynomial being assigned the fields of <i>f</i>.
     */
    inline void swap( BinaryPolynomial & f , BinaryPolynomial & g ) {

        BinaryPolynomial::swap(f,g);
    }

    /**
     * @brief
     *           Computes the sum of two binary polynomials.
     *
     * @details
     *           The method replaces the polynomial <i>f</i>
     *           by the sum of the two polynomials <i>g</i>
     *           and <i>h</i>.
     *
     * @param f
     *           Sum of <i>g</i> and <i>h</i>.
     *
     * @param g
     *           First summand.
     *
     * @param h
     *           Second summand.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void add
    ( BinaryPolynomial & f ,
      const BinaryPolynomial & g ,
      const BinaryPolynomial & h ) {

        BinaryPolynomial::add(f,g,h);
    }

    /**
     * @brief
     *           Computes the difference of two binary polynomials.
     *
     * @details
     *           The method replaces the polynomial <i>f</i>
     *           by the difference of the two polynomials <i>g</i>
     *           and <i>h</i>.
     *
     *           Note that, since we are computing with binary
     *           arithmetic, computing the difference between
     *           two polynomials is equivalent to addition.
     *
     * @param f
     *           Difference of <i>g</i> and <i>h</i>.
     *
     * @param g
     *           Minuend.
     *
     * @param h
     *           Subtrahend.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void sub
    ( BinaryPolynomial & f ,
      const BinaryPolynomial & g ,
      const BinaryPolynomial & h ) {

        BinaryPolynomial::sub(f,g,h);
    }

    /**
     * @brief
     *           Computes the produce of two binary polynomials.
     *
     * @details
     *           The method replaces the polynomial <i>f</i>
     *           by the product of the two polynomials <i>g</i>
     *           and <i>h</i>.
     *
     * @param f
     *           Product of <i>g</i> and <i>h</i>.
     *
     * @param g
     *           First factor.
     *
     * @param h
     *           Second factor.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void mul
    ( BinaryPolynomial & f ,
      const BinaryPolynomial & g ,
      const BinaryPolynomial & h ) {

        BinaryPolynomial::mul(f,g,h);
    }

    /**
     * @brief
     *           Computes a Euclidean division with remainder.
     *
     * @details
     *           Given two polynomials <i>a</i> and a non-zero <i>b</i>
     *           there exits unique polynomials <i>q</i> and <i>r</i>
     *           such that
     *           \f$
     *            a=q\cdot b+r
     *           \f$
     *           with \f$\deg(r)<\deg(b)\f$. These two polynomials
     *           are determined by the method.
     *
     * @param q
     *           The quotient of the Euclidean division.
     *
     * @param r
     *           The remainder of the Euclidean division.
     *
     * @param a
     *           Numerator.
     *
     * @param b
     *           Denominator.
     *
     * @warning
     *           If <i>b</i> is zero, an error message is printed
     *           to <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If <i>q</i> and <i>r</i> are the same references,
     *           an error message is printed to <code>stderr</code>
     *           and the program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void divRem
    ( BinaryPolynomial & q , BinaryPolynomial & r ,
      const BinaryPolynomial & a , const BinaryPolynomial & b ) {

        BinaryPolynomial::divRem(q,r,a,b);
    }

    /**
     * @brief
     *           Computes the inverse of a polynomial modulo
     *           another polynomial.
     *
     * @details
     *           The method attempts to find a polynomial <i>s</i>
     *           of degree smaller than \f$\deg(b)\f$ such that
     *           \f[
     *            1=c(X)\cdot a(X)~rem~b.
     *           \f]
     *           If such a polynomial exists, the function
     *           returns <code>true</code> and replaces the content
     *           of <i>c</i> correspondingly; otherwise, if the inverse
     *           does not exist, the function returns <code>false</code>
     *           and leaves the content of <i>c</i> unchanged.
     *
     * @param c
     *           Output of the inverse of <i>a</i> modulo <i>b</i>.
     *
     * @param a
     *           Input polynomial of which the inverse is computed.
     *
     * @param b
     *           Input modulus polynomial.
     *
     * @return
     *           <code>true</code> if the inverse have been determined
     *           successfully, in particular, if it exists; otherwise
     *           the function returns <code>false</code>.
     *
     * @warning
     *           If <i>a</i> and <i>b</i> are both zero, the program
     *           prints an error message to <code>stderr</code> and
     *           exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline bool invMod
    ( BinaryPolynomial & c ,
      const BinaryPolynomial & a , const BinaryPolynomial & b ) {

        return BinaryPolynomial::invMod(c,a,b);
    }

    /**
     * @brief
     *           Computes the quotient of two polynomials via Euclidean
     *           division.
     *
     * @details
     *           The method calls the \link divRem()\endlink
     *           method and dismisses the remainder.
     *
     * @param q
     *           Quotient.
     *
     * @param a
     *           Numerator.
     *
     * @param b
     *           Denominator.
     *
     * @see divRem()
     *
     * @warning
     *           If <i>b</i> is zero, an error message is printed
     *           to <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void div
    ( BinaryPolynomial & q , const BinaryPolynomial & a ,
      const BinaryPolynomial & b ) {

        BinaryPolynomial::div(q,a,b);
    }

    /**
     * @brief
     *           Computes the remainder of two polynomials via Euclidean
     *           division.
     *
     * @details
     *           The method calls the \link divRem()\endlink
     *           method and dismisses the quotient.
     *
     * @param r
     *           Remainder.
     *
     * @param a
     *           Numerator.
     *
     * @param b
     *           Denominator.
     *
     * @see divRem()
     *
     * @warning
     *           If <i>b</i> is zero, an error message is printed
     *           to <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void rem
    ( BinaryPolynomial & r , const BinaryPolynomial & a ,
      const BinaryPolynomial & b ) {

        BinaryPolynomial::rem(r,a,b);
    }

    /**
     * @brief
     *          Computes the greatest common divisor of two binary
     *          polynomials.
     *
     * @param g
     *          The greatest common divisors of <i>a</i> and <i>b</i>.
     *
     * @param a
     *          First polynomial.
     *
     * @param b
     *          Second polynomial.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void gcd
    ( BinaryPolynomial & g ,
      const BinaryPolynomial & a , const BinaryPolynomial & b ) {

        BinaryPolynomial::gcd(g,a,b);
    }


    /**
     * @brief
     *           Performs the extended Euclidean algorithm.
     *
     * @details
     *           The method determines the greatest common divisor
     *           <i>g</i> of the two polynomials <i>a</i> and <i>b</i>
     *           and co-prime polynomials <i>s</i> and <i>t</i> such that
     *           \f[
     *            g=s\cdot a+t\cdot b.
     *           \f]
     *
     * @param g
     *           Output of the greatest common divisor of <i>a</i>
     *           and <i>b</i>.
     *
     * @param s
     *           Output polynomial as described above.
     *
     * @param t
     *           Output polynomial as described above.
     *
     * @param a
     *           First input polynomial.
     *
     * @param b
     *           Second input polynomial.
     *
     * @warning
     *           If <i>g</i>, <i>s</i> or <i>t</i> are of the same
     *           reference, an error message is printed to
     *           <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If <i>a</i> and <i>b</i> are both zero, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void xgcd
    ( BinaryPolynomial & g ,
      BinaryPolynomial & s , BinaryPolynomial & t ,
      const BinaryPolynomial & a , const BinaryPolynomial & b ) {

        BinaryPolynomial::xgcd(g,s,t,a,b);
    }

    /**
     * @brief
     *            Computes the partial greatest common divisor of two
     *            polynomials.
     *
     * @details
     *            This methods performs the traditional extended euclidean
     *            algorithm, i.e., determines \f$s,t\f$ such that
     *            \f$g=s\cdot a+t\cdot b\f$, until the degree of
     *            \f$g\f$ becomes smaller than \f$d\f$.
     *
     * @param g
     *            Will contain the partial greatest common divisor
     *            of <i>a</i> and <i>b</i>.
     *
     * @param s
     *            see details
     * @param t
     *            see details
     * @param a
     *            see details
     * @param b
     *            see details
     * @param d
     *            see details
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void pgcd
    ( BinaryPolynomial & g ,
      BinaryPolynomial & s ,
      BinaryPolynomial & t ,
      const BinaryPolynomial & a ,
      const BinaryPolynomial & b ,
      int d ) {

        BinaryPolynomial::pgcd(g,s,t,a,b,d);
    }

    /**
     * @brief
     *            Computes the <i>(n,n-k)</i>-Pade approximant to the
     *            polynomial.
     *
     * @details
     *            The method determines co-prime polynomials <i>r</i> and
     *            <i>t</i> with \f$\deg(r)<k\f$ and \f$\deg(t)\leq n-k\f$
     *            such that
     *            \f[
     *             \frac{r(X)}{t(X)}=g(X)~mod~X^n
     *            \f]
     *            where \f$t(X)\f$ is not divisible by \f$X\f$.
     *
     *            For more details on <em>Pade approximation</em> we refer
     *            to Section 5.9 in
     *            <ul>
     *             <li><b>[vzGth]</b> v.z. Gathen and Gerhard (2003).
     *              <i>Modern Computer Algebra</i>. Cambridge University
     *              Press, Cambridge (UK), 2nd edition.
     *             </li>
     *            </ul>
     *
     * @param r
     *            Output numerator polynomial.
     *
     * @param t
     *            Output denominator polynomial.
     *
     * @param g
     *            Input polynomial.
     *
     * @param n
     *            Bounds the degree of <i>r</i>.
     *
     * @param k
     *            Bounds the degree of <i>t</i>, i.e., such that
     *            \f$\deg(t)\leq n-k\f$.
     *
     * @warning
     *            If <i>n</i> or <i>k</i> is negative, the degree
     *            of <i>g</i> is greater than or equals <i>n</i>, or
     *            if <i>k>n</i>, an error message is written to
     *            <code>stderr</code> and the program exits with
     *            status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void pade
    ( BinaryPolynomial & r , BinaryPolynomial & t ,
      const BinaryPolynomial & g , int n , int k ) {

        BinaryPolynomial::pade(r,t,g,n,k);
    }

    /**
     * @brief
     *           Computes the minimal polynomial of a polynomial modulo
     *           another polynomial.
     *
     * @details
     *           The minimal polynomial of <i>g</i> modulo <i>f</i> is
     *           defined as the non-zero polynomial <i>h</i> of minimal
     *           degree such that
     *           \f[
     *            h(g(X))~mod~f(X)=0.
     *           \f]
     *
     * @param h
     *           Output to contain the minimal polynomial of <i>g</i>
     *           modulo <i>f</i>.
     *
     * @param g
     *           Input polynomial.
     *
     * @param f
     *           Input modulus polynomial; must be of degree greater
     *           than 1.
     *
     * @warning
     *           If <i>f</i> is of degree smaller than 1, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void minPolyMod
    ( BinaryPolynomial & h ,
      const BinaryPolynomial & g , const BinaryPolynomial & f ) {

        BinaryPolynomial::minPolyMod(h,g,f);
    }

    /**
     * @brief
     *          Computes the composition of two binary polynomials.
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
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void eval
    ( BinaryPolynomial & h ,
      const BinaryPolynomial & f ,
      const BinaryPolynomial & g ) {

        BinaryPolynomial::eval(h,f,g);
    }

    /**
     * @brief
     *           Computes the composition of two polynomials being
     *           reduced by a modulus.
     *
     * @details
     *           The method performs <i>h(X)=f(g(X))</i> mod <i>m</i>.
     *
     * @param h
     *           Modular composition.
     *
     * @param f
     *           The outer polynomial.
     *
     * @param g
     *           The inner polynomial.
     *
     * @param m
     *           Modulus.
     *
     * @warning
     *           If <i>m</i> is zero, an error message is printed
     *           to <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     *
     */
    inline void evalMod
    ( BinaryPolynomial & h , const BinaryPolynomial & f ,
      const BinaryPolynomial & g , const BinaryPolynomial & m ) {

        BinaryPolynomial::evalMod(h,f,g,m);
    }

    /**
     * @brief
     *          Computes the square root of a polynomial.
     *
     * @details
     *          If the input polynomial <i>f</i> is the square
     *          of another polynomial <i>g</i>, i.e., \f$f=g^2\f$, then
     *          the method will correctly determined <i>g</i>.
     *
     * @param g
     *          On output square root of the input polynomial <i>f</i>.
     *
     * @param f
     *          On input the square of the polynomial <i>g</i>.
     *
     * @see isSquare()
     *
     * @warning
     *          If <i>f</i> is not the square of a polynomial, the
     *          result of the method is undocumented.
     *
     * @warning
     *          If not enough memory could be provided, an error
     *          message is printed to <code>stderr</code> and the
     *          program exits with status 'EXIT_FAILURE'.
     */
    inline void squareRoot
        ( BinaryPolynomial & g , const BinaryPolynomial & f ) {

        BinaryPolynomial::squareRoot(g,f);
    }

    /**
     * @brief
     *            Prints a text representation of a polynomial to the
     *            specified output stream.
     *
     * @param out
     *            The output stream.
     *
     * @param f
     *            The polynomial of which a text representation is
     *            printed to <code>out</code>.
     *
     * @return
     *            A reference to <code>out</code> after the text
     *            representation has been written.
     */
    THIMBLE_DLL std::ostream & operator<<
    ( std::ostream & out , const BinaryPolynomial & f );
}

#endif /* THIMBLE_BINARYPOLYNOMIAL_H */
