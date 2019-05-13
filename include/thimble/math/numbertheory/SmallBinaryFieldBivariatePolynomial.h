/*
 *  THIMBLE --- Research Library for Development and Analysis of
 *  Fingerprint-Based Biometric Cryptosystems.
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
 * @file SmallBinaryFieldBivariatePolynomial.h
 *
 * @brief
 *            Provides a class of which objects represent bivariate
 *            polynomials over small binary galois fields.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_SMALLBINARYFIELDBIVARIATEPOLYNOMIAL_H_
#define THIMBLE_SMALLBINARYFIELDBIVARIATEPOLYNOMIAL_H_

#include <new>

#include <thimble/dllcompat.h>
#include <thimble/math/numbertheory/SmallBinaryField.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Objects of this class represent bivariate polynomials
	 *            with coefficients in a small binary galois field.
	 */
	class THIMBLE_DLL SmallBinaryFieldBivariatePolynomial {

	private:

		/**
		 * @brief
		 *            Pointer to the galois field from which the coefficients
		 *            of this bivariate polynomial are chosen.
		 */
		const SmallBinaryField *gfPtr;

		/**
		 * @brief
		 *            Y-coefficients of this polynomial.
		 *
		 * @details
		 *            Bivariate polynomials are polynomials in two variables
		 *            having coefficients in a field <i>F</i>, say, and are
		 *            denoted by the set
		 *            \f[
		 *             F[X,Y]=\{\sum_{i=0}^{d_x}\sum_{j=0}^{d_y}c_{i,j}X^iY^j\}.
		 *            \f]
		 *            Internally, a %SmallBinaryFieldBivariatePolynomial
		 *            object does use the representation \f$F[X][Y]\f$,
		 *            i.e.,
		 *            \f[
		 *             F[X][Y] = \{\sum_{j=0}^{d_y}c_j(X)\cdot Y^j\}
		 *            \f]
		 *            where the \f$c_j(X)\f$ are polynomials in \f$F[X]\f$
		 *            such that
		 *            \f[
		 *             c_j(X)=\sum_{i=0}^{d_x}c_i\cdot X^i.
		 *            \f]
		 *
		 * @see getCoeffY()
		 * @see setCoeffY()
		 */
		SmallBinaryFieldPolynomial *coefficientsY;

		/**
		 * @brief
		 *            The Y-degree of this bivariate polynomial.
		 *
		 * @details
		 *            Let
		 *            \f[
		 *          \sum_{j=0}^{d_y}c_j(X)\cdot Y^j
		 *            \f]
		 *            be the polynomial represented by this object
		 *            where \f$c_j(X)\f$ are the polynomials successively
		 *            stored in \link coefficientsY\endlink such that
		 *            \f$c_{d_y}(X)\neq 0\f$. Then <code>degreeY=</code>
		 *            \f$d_y\f$.
		 *
		 * @see degY()
		 */
		int degreeY;

		/**
		 * @brief
		 *            The number of uni-variate polynomials (in X)
		 *            that \link coefficientsY\endlink is able
		 *            to store without requiring reallocation.
		 *
		 * @attention
		 *            <b>Note:</b> the capacity can be larger
		 *            than \link degreeY\endlink+1.
		 */
		int capacity;

		/**
		 * @brief
		 *            Ensures that \link degreeY\endlink correctly encodes
		 *            the Y-degree of this polynomial.
		 */
		void normalize();

		/**
		 * @brief
		 *            Computes the product between two bivariate polynomials
		 *            assuming that the specified references are different
		 *            and the bivariate polynomial that they represents have
		 *            been specified with the same coefficient field.
		 *
		 * @param h
		 *            Product of <code>f</code> and <code>g</code>.
		 *
		 * @param f
		 *            First factor.
		 *
		 * @param g
		 *            Second factor.
		 *
		 * @warning
		 *            If <code>h</code>, <code>f</code> or <code>g</code> have
		 *            have been different over distinct finite fields or if
		 *            they are not of pairwise different reference, the program
		 *            runs into undocumented behavior.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		static void mulUncheck
		( SmallBinaryFieldBivariatePolynomial & h ,
		  const SmallBinaryFieldBivariatePolynomial & f ,
		  const SmallBinaryFieldBivariatePolynomial & g );

	public:

		/**
		 * @brief
		 *            Ensures that this bivariate polynomial can hold at least
		 *            the specified number of Y-coefficients without requiring
		 *            reallocation.
		 *
		 * @details
		 *            If during computations a maximum number of Y-coefficients
		 *            can be predicted, this maximum can be specified via
		 *            this method to safe time-consuming internal reallocations.
		 *
		 *            Note that, if this object can already hold at least the
		 *            specified number of Y-coefficients, calling this method
		 *            is without effect. In particular, the specified capacity
		 *            is allowed to be zero or negative.
		 *
		 * @param capacity
		 *            The number of Y-coefficients this bivariate polynomial
		 *            will be able to hold.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		void ensureCapacity( int capacity );

		/**
		 * @brief
		 *            Assignment operator (procedural version).
		 *
		 * @param f
		 *            The bivariate polynomial of which this object is assigned
		 *            with a copy.
		 *
		 * @warning
		 *            If this object's coefficient field is different from
		 *            the field of <code>f</code>, then the program may
		 *            print an error message to <code>stderr</code> and
		 *            exit with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		void assign( const SmallBinaryFieldBivariatePolynomial & f );

		/**
		 * @brief
		 *            Essentially, a copy constructor of an uni-variate
		 *            polynomial in X represented as a bivariate polynomial
		 *            being constant in Y.
		 *
		 * @details
		 *            If <code>c0</code> represents the univariate polynomial
		 *            \f[
		 *             c_0(X)=\sum_{i=0}^{d_x}c_{i,0}X^i,
		 *            \f]
		 *            then the bivariate polynomial
		 *            \f[
		 *             f(X,Y)=c_0(X)
		 *            \f]
		 *            that is constant in Y is created.
		 *
		 * @param c0
		 *            The univariate polynomial of which a copy as a bivariate
		 *            polynomial is created.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
        SmallBinaryFieldBivariatePolynomial
        ( const SmallBinaryFieldPolynomial & c0 );

        /**
         * @brief
         *            Destructor.
         *
         * @details
         *            Frees all the memory that has been allocated to
         *            represent this bivariate polynomial.
         */
        ~SmallBinaryFieldBivariatePolynomial();

        /**
         * @brief
         *            Creates a zero bivariate polynomial of which
         *            coefficients are zero.
         *
         * @param gf
         *            The finite field from which the (initially zero)
         *            coefficients are drawn.
         */
		inline SmallBinaryFieldBivariatePolynomial
		( const SmallBinaryField & gf ) {
			this->gfPtr = &gf;
			this->coefficientsY = NULL;
			this->degreeY = -1;
			this->capacity = 0;
		}

		/**
		 * @brief
		 *            Copy constructor.
		 *
		 * @details
		 *            Constructs a bivariate polynomial with coefficients
		 *            the same finite field as <code>f</code> that represents
		 *            a copy of <code>f</code>
		 *
		 * @param f
		 *            The bivariate polynomial of which the copy is created.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		inline SmallBinaryFieldBivariatePolynomial
		( const SmallBinaryFieldBivariatePolynomial & f ) {
			this->gfPtr = f.gfPtr;
			this->coefficientsY = NULL;
			this->degreeY = -1;
			this->capacity = 0;
			assign(f);
		}

		/**
		 * @brief
		 *            Assignment operator.
		 *
		 * @details
		 *            The overloaded '='-operator is implemented by wrapping
		 *            around the \link assign()\endlink method.
		 *
		 * @param f
		 *            The bivariate polynomial of which this object is assigned
		 *            with a copy.
		 *
		 * @return
		 *            A reference to this object (after assignment)
		 *
		 * @see assign()
		 *
		 * @warning
		 *            If this object's coefficient field is different from
		 *            the field of <code>f</code>, then the program may
		 *            print an error message to <code>stderr</code> and
		 *            exit with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		inline SmallBinaryFieldBivariatePolynomial& operator=
		( const SmallBinaryFieldBivariatePolynomial & f ) {
			assign(f);
			return *this;
		}

		/**
		 * @brief
		 *            Swaps the representation (content) of this polynomial
		 *            with the representation (content) of the specified
		 *            polynomial.
		 *
		 * @details
		 *            After the method has been called this polynomial's
		 *            representation will be exchanged with the representation
		 *            of the polynomial specified by <code>f</code>
		 *            (including the references the underlying finite fields,
		 *             if different).
		 *
		 * @param f
		 *            The bivariate polynomial of which representation is
		 *            exchanged with the representation of this bivariate
		 *            polynomial.
		 */
		inline void swap( SmallBinaryFieldBivariatePolynomial & f ) {

			const SmallBinaryField *gfPtr;
			SmallBinaryFieldPolynomial *coefficientsY;
			int degreeY;
			int capacity;

			gfPtr         = this->gfPtr;
			coefficientsY = this->coefficientsY;
			degreeY       = this->degreeY;
			capacity      = this->capacity;

			this->gfPtr         = f.gfPtr;
			this->coefficientsY = f.coefficientsY;
			this->degreeY       = f.degreeY;
			this->capacity      = f.capacity;

			f.gfPtr         = gfPtr;
			f.coefficientsY = coefficientsY;
			f.degreeY       = degreeY;
			f.capacity      = capacity;
		}

		/**
		 * @brief
		 *            Accesses the finite field from which the coefficients
		 *            of this bivariate polynomial are drawn.
		 *
		 * @return
		 *            A constant reference to this polynomial's underlying
		 *            finite field.
		 */
		inline const SmallBinaryField & getField() const {
			return this->gfPtr[0];
		}

		/**
		 * @brief
		 *            Accesses the Y-degree of this bivariate polynomial.
		 *
		 * @details
		 *            Write
		 *            \f[
		 *             f(X,Y)=\sum_{j=0}^{d_y}c_j(X)\cdot Y^{d_y}
		 *            \f]
		 *            where \f$c_{d_y}\neq 0\f$. Then \f$d_y\f$ is defined
		 *            to be the Y-degree of the bivariate polynomial \f$f\f$.
		 *
		 * @return
		 *            The Y-degree of this bivariate polynomial
		 */
		inline int degY() const {
			return this->degreeY;
		}

		/**
		 * @brief
		 *            Returns the X-degree of the <i>j</i>th Y-coefficient.
		 *
		 * @details
		 *            Write the polynomial as
		 *            \f[
		 *             \sum_{j}c_j(X)\cdot Y^j
		 *            \f]
		 *            where the \f$c_j(X)\f$ are polynomials in the
		 *            indeterminate \f$X\f$ with coefficients in the finite
		 *            field. The result is the X-degree of \f$c_j(X)\f$.
		 *
		 * @see SmallBinaryFieldPolynomial::deg()
		 *
		 * @return
		 *            The X-degree of the <i>j</i>th Y-coefficient.
		 */
		inline int degX( int j ) const {
			if ( j > degY() ) {
				return -1;
			}
			return this->coefficientsY[j].deg();
		}

		/**
		 * @brief
		 *            Tests if this bivariate polynomial is zero.
		 *
		 * @return
		 *            <code>true</code> if this polynomial is identical
		 *            to zero or, otherwise, <code>false</code>.
		 */
		inline bool isZero() const {
			return this->degreeY < 0;
		}

		/**
		 * @brief
		 *            Computes the <i>(a,b)</i>-weighted degree of this
		 *            bivariate polynomial.
		 *
		 * @details
		 *            The <i>(a,b)</i>-weighted degree of a bivariate
		 *            polynomial
		 *            \f[
		 *             f(X,Y) = \sum_{i,j}c_{i,j}X^iY^j
		 *            \f]
		 *            is defined as the expression
		 *            \f[
		 *             \max_{c_{i,j}\neq 0}\{a\cdot i+b\cdot j\}.
		 *            \f]
		 *
		 * @return
		 *            The <i>(a,b)</i>-weighted degree of this bivariate
		 *            polynomial.
		 */
		std::pair<int,int> deg( int a , int b ) const;

		/**
		 * @brief
		 *            Access the <i>i</i>th X and the <i>j</i>th
		 *            Y-coefficient.
		 *
		 * @details
		 *            Write
		 *            \f[
		 *             f(X,Y) = \sum_{i,j}c_{i,j}X^iY^j
		 *            \f]
		 *            where \f$c_{i,j}\f$ denote the coefficient
		 *            in the finite field of this bivariate polynomial.
		 *            Then this function returns the representation
		 *            for \f$c_{i,j}\f$.
		 *
		 *            Note, this function returns 0 if <code>i</code> or
		 *            <code>j</code> is negative or if they exceed the
		 *            X or Y-degree, respectively.
		 *
		 * @param i
		 *            The index of the X-coefficient.
		 *
		 * @param j
		 *            The index of the Y-coefficient.
		 *
		 * @return
		 *            The <i>i</i>th X and the <i>j</i>th Y-coefficient.
		 */
		uint32_t getCoeff( int i , int j ) const;

		/**
		 * @brief
		 *            Access the <i>j</i>th Y-coefficient of this
		 *            bivariate polynomial.
		 *
		 * @details
		 *            Write
		 *            \f[
		 *             f(X,Y) = \sum_{j}c_{j}(X)Y^j
		 *            \f]
		 *            for this polynomial where the \f$c_j(X)\f$ denote
		 *            the univariate Y-coefficient polynomials in X. This
		 *            function returns \f$c_j(X)\f$.
		 *
		 * @param j
		 *            The index of Y-coefficient.
		 *
		 * @return
		 *            The <i>j</i>th Y-coefficient of this bivariate
		 *            polynomial.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		inline SmallBinaryFieldPolynomial getCoeffY( int j ) const {
			if ( j < 0 || j > degY() ) {
				return SmallBinaryFieldPolynomial(*gfPtr);
			}
			return this->coefficientsY[j];
		}

		/**
		 * @brief
		 *           Determine whether this bivariate polynomial is equals to
		 *           the given bivariate polynomial.
		 *
		 * @details
		 *           If this polynomial's coefficients are the same
		 *           as the coefficients of <code>f</code> the function
		 *           return <code>true</code> and <code>false</code>
		 *           otherwise.
		 *
		 * @return
		 *           <code>true</code> if this polynomial's
		 *           coefficients are the same as the coefficients of
		 *           <code>f</code> and <code>false</code> otherwise.
		 *
		 * @warning
		 *           If this polynomial is defined over a finite field which
		 *           is different to the finite field of <code>f</code>
		 *           (i.e. if the corresponding reference is different) then
		 *           an error message is printed to <code>stderr</code>
		 *           and the program exits with status 'EXIT_FAILURE'.
		 */
		bool equals( const SmallBinaryFieldBivariatePolynomial & f) const;

		/**
		 * @brief
		 *           Makes this bivariate polynomial the constant zero
		 *           polynomial.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		inline void setZero() {
			ensureCapacity(1);
			this->coefficientsY[0].setZero();
			this->degreeY = -1;
		}

		/**
		 * @brief
		 *            Sets this polynomial equals \f$Y\f$.
		 *
		 * @details
		 *            This polynomial will be set to
		 *            \f[
		 *             f(X,Y)\leftarrow Y.
		 *            \f]
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		inline void setY() {
			ensureCapacity(2);
			this->coefficientsY[0].setZero();
			this->coefficientsY[1].setOne();
			this->degreeY = 1;
		}

		/**
		 * @brief
		 *            Replaces the <i>i</i>th Y-coefficient.
		 *
		 * @details
		 *            Write
		 *            \f[
		 *             \sum_{j}c_j(X)Y^j
		 *            \f]
		 *            for this bivariate polynomial. This methods
		 *            performs
		 *            \f[
		 *             c_i(X)\leftarrow c(X)
		 *            \f]
		 *            where \f$c(X)\f$ is specified by <code>c</code>.
		 *
		 * @param i
		 *            The position at which the Y-coefficient is updated.
		 *
		 * @param c
		 *            The new <i>i</i>th Y-coefficient.
		 *
		 * @warning
		 *            If <i>i</i> is negative, then an error message will
		 *            be printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <code>c</code> is defined over a different finite
		 *            field than this bivariate polynomial, then an error
		 *            message will be printed to <code>stderr</code> and
		 *            the program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		void setCoeffY( int i , const SmallBinaryFieldPolynomial & c );

		/**
		 * @brief
		 *            Updates the <i>i</i>th X and <i>j</i>th Y coefficient.
		 *
		 * @param i
		 *            The index of the X coefficient.
		 *
		 * @param j
		 *            The index of the Y coefficient.
		 *
		 * @param c
		 *            The new coefficient.
		 *
		 * @warning
		 *            If <i>i</i> or <i>j</i> negative, then an error message
		 *            will be printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <code>c</code> does not represent a valid element in
		 *            the finite field associated with this bivariate
		 *            polynomial, then an error message will be printed to
		 *            <code>stderr</code> and the program exits with status
		 *            'EXIT_FAILURE'
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		void setCoeff( int i , int j , uint32_t c );

		/**
		 * @brief
		 *            Evaluates this polynomial at specified position.
		 *
		 * @details
		 *            If this polynomial is denoted by \f$f(X,Y)\f$ in the
		 *            indeterminate \f$X\f$ and \f$Y\f$, then this function
		 *            returns \f$f(x,y)\f$.
		 *
		 * @param x
		 *            First position component.
		 *
		 * @param y
		 *            Second position component.
		 *
		 * @return
		 *            This bivariate polynomial evaluated at <i>(x,y)</i>.
		 *
		 * @warning
		 *            If <code>x</code> or <code>y</code> do not represent
		 *            valid elements in the finite field associated with this
		 *            bivariate polynomial, then this function may print an
		 *            error message to <code>stderr</code> and exits with
		 *            status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		uint32_t eval( uint32_t x , uint32_t y ) const;

		/**
		 * @brief
		 *            Computes the representation of a bivariate polynomial
		 *            after exchanging its two indeterminate symbols.
		 *
		 * @param g
		 *            The input polynomial \f$g(X,Y)\f$.
		 *
		 * @param f
		 *            The output polynomial which will be
		 *            equals \f$f(X,Y)=g(Y,X)\f$.
		 *
		 * @see \link thimble::swapXY() swapXY()\endlink
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		static void swapXY
		( SmallBinaryFieldBivariatePolynomial & f ,
		  const SmallBinaryFieldBivariatePolynomial & g );

		/**
		 * @brief
		 *            Determines the leading term of the highest
		 *            <i>(a,b)</i>-weighted coffient.
		 *
		 * @details
		 *            Write for this polynomial
		 *            \f[
		 *             \sum_{i,j}c_{i,j}X^iY^j.
		 *            \f]
		 *            This function determines <i>(i,j)</i> for which the
		 *            expression <i>a*i+b*j</i> is maximal and then sets
		 *            <code>lt</code> the a representation of
		 *            \f[
		 *             c_{i,j}X^iY^j.
		 *            \f]
		 *            The pair <i>(i,j)</i> will be furthermore returned.
		 *
		 * @param lt
		 *            Reference to which the leading term, according to the
		 *            <i>(a,b)</i>-weighted degree, will be stored.
		 *
		 * @param f
		 *            Input bivariate polynomial from which the leading
		 *            term is computed.
		 *
		 * @param a
		 *            First weight.
		 *
		 * @param b
		 *            Second weight.
		 *
		 * @see \link thimble::leadTerm(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,int,int) leadTerm()\endlink
		 * @see deg()
		 *
		 * @return
		 *            The indices <i>(i,j)</i> at which this polynomial's
		 *            cofficient is non-zero and for which the expression
		 *            <i>a*i+b*j</i> is maximal.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		static std::pair<int,int> leadTerm
		( SmallBinaryFieldBivariatePolynomial & lt ,
		  const SmallBinaryFieldBivariatePolynomial & f ,
		  int a , int b );

		/**
		 * @brief
		 *            Essentially, multiplies a bivariate polynomial with
		 *            the factor \f$Y^n\f$.
		 *
		 * @details
		 *            Write
		 *            \f[
		 *             g(X,Y) = \sum_{i,j\geq 0}c_{i,j}X^iY^j.
		 *            \f]
		 *            This method computes
		 *            \f[
		 *             f(X,Y) = Y^n\cdot g(X,Y)
		 *            \f]
		 *            if \f$n\f$ is positive and
		 *            \f[
		 *             f(X,Y) = \sum_{i\geq 0,j\geq -n}c_{i,j}X^iY^{j+n}
		 *            \f]
		 *            if \f$n\f$ is negative.
		 *
		 * @param f
		 *            The output polynomial being the product.
		 *
		 * @param g
		 *            The input factor polynomial.
		 *
		 * @param n
		 *            The exponent of the factor monomial.
		 *
		 * @see \link thimble::leftShiftY(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,int) leftShiftY()\endlink
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		static void leftShiftY
		( SmallBinaryFieldBivariatePolynomial & f ,
		  const SmallBinaryFieldBivariatePolynomial & g ,
		  int n );

		/**
		 * @brief
		 *            Essentially, computes the quotient by a division
		 *            of a bivariate polynomial divided by \f$Y^n\f$.
		 *
		 * @details
		 *            Write
		 *            \f[
		 *             g(X,Y) = \sum_{i,j\geq 0}c_{i,j}X^iY^j.
		 *            \f]
		 *            This method computes the polynomial
		 *            \f[
		 *             f(X,Y) = \sum_{i\geq 0,j\geq n}c_{i,j}X^iY^{j-n}.
		 *            \f]
		 *            for positive \f$n\f$. If \f$n\f$ is negative, then
		 *            \f[
		 *             f(X,Y) = Y^{-n}\cdot g(X,Y)
		 *            \f]
		 *            is computed.
		 *
		 * @param f
		 *            Output bivariate polynomial the will contain the
		 *            quotient.
		 *
		 * @param g
		 *            Input polynomial divided by \f$Y^n\f$.
		 *
		 * @param n
		 *            The exponent of the divisor monomial.
		 *
		 * @see \link thimble::rightShiftY(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,int) rightShiftY()\endlink
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		static void rightShiftY
		( SmallBinaryFieldBivariatePolynomial & f ,
		  const SmallBinaryFieldBivariatePolynomial & g ,
		  int n );

		/**
		 * @brief
		 *            Essentially, multiplies a bivariate polynomial with
		 *            the factor \f$X^n\f$.
		 *
		 * @details
		 *            Write
		 *            \f[
		 *             g(X,Y) = \sum_{i,j\geq 0}c_{i,j}X^iY^j.
		 *            \f]
		 *            This method computes
		 *            \f[
		 *             f(X,Y) = X^n\cdot g(X,Y)
		 *            \f]
		 *            if \f$n\f$ is positive and
		 *            \f[
		 *             f(X,Y) = \sum_{i\geq -n,j\geq 0}c_{i,j}X^{i+n}Y^{j}
		 *            \f]
		 *            if \f$n\f$ is negative.
		 *
		 * @param f
		 *            The output polynomial being the product.
		 *
		 * @param g
		 *            The input factor polynomial.
		 *
		 * @param n
		 *            The exponent of the factor monomial.
		 *
		 * @see \link thimble::leftShiftX(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,int) leftShiftX()\endlink
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		static void leftShiftX
		( SmallBinaryFieldBivariatePolynomial & f ,
		  const SmallBinaryFieldBivariatePolynomial & g ,
		  int n );

		/**
		 * @brief
		 *            Essentially, computes the quotient by a division
		 *            of a bivariate polynomial divided by \f$X^n\f$.
		 *
		 * @details
		 *            Write
		 *            \f[
		 *             g(X,Y) = \sum_{i,j\geq 0}c_{i,j}X^iY^j.
		 *            \f]
		 *            This method computes the polynomial
		 *            \f[
		 *             f(X,Y) = \sum_{i\geq n,j\geq 0}c_{i,j}X^{i-n}Y^{j}.
		 *            \f]
		 *            for positive \f$n\f$. If \f$n\f$ is negative, then
		 *            \f[
		 *             f(X,Y) = X^{-n}\cdot g(X,Y)
		 *            \f]
		 *            is computed.
		 *
		 * @param f
		 *            Output bivariate polynomial the will contain the
		 *            quotient.
		 *
		 * @param g
		 *            Input polynomial divided by \f$X^n\f$.
		 *
		 * @param n
		 *            The exponent of the divisor monomial.
		 *
		 * @see \link thimble::rightShiftX(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,int) rightShiftX()\endlink
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		static void rightShiftX
		( SmallBinaryFieldBivariatePolynomial & f ,
		  const SmallBinaryFieldBivariatePolynomial & g ,
		  int n );

		/**
		 * @brief
		 *            Computes the sum of two bivariate polynomials.
		 *
		 * @details
		 *            This method computes
		 *            \f[
		 *             h(X,Y)=f(X,Y)+g(X,Y).
		 *            \f]
		 *
		 * @param h
		 *            Output polynomial that will contain the sum.
		 *
		 * @param f
		 *            First summand.
		 *
		 * @param g
		 *            Second summand.
		 *
		 * @see \link thimble::add(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&) add()\endlink
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		static void add
		( SmallBinaryFieldBivariatePolynomial & h ,
		  const SmallBinaryFieldBivariatePolynomial & f ,
		  const SmallBinaryFieldBivariatePolynomial & g );

		/**
		 * @brief
		 *            Adds an univariate polynomial to the lowest
		 *            Y-coefficient of a bivariate polynomial.
		 *
		 * @details
		 *            This method computes
		 *            \f[
		 *             h(X,Y)=f(X,Y)+s(X).
		 *            \f]
		 *
		 * @param h
		 *            Output polynomial that will contain the sum.
		 *
		 * @param f
		 *            First bivariate summand.
		 *
		 * @param s   Second univariate summand.
		 *
		 * @see \link thimble::add(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldPolynomial&) add()\endlink
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		static void add
		( SmallBinaryFieldBivariatePolynomial & h ,
		  const SmallBinaryFieldBivariatePolynomial & f ,
		  const SmallBinaryFieldPolynomial & s );

		/**
		 * @brief
		 *            Computes difference between two bivariate polynomials
		 *            which, actually, is equivalent to the computation of
		 *            their sum in binary fields.
		 *
		 * @details
		 *            This method computes
		 *            \f[
		 *             h(X,Y)=f(X,Y)-g(X,Y)
		 *            \f]
		 *            which is equivalent to the computation of their sum.
		 *            Nonetheless, we provide subtraction and negation
		 *            methods explicitly in order to allow the user of
		 *            this library to produce possibly better readable
		 *            code.
		 *
		 * @param h
		 *            Output polynomial that will contain the sum.
		 *
		 * @param f
		 *            Minuend.
		 *
		 * @param g
		 *            Subtrahend.
		 *
		 * @see \link thimble::sub(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&) add()\endlink
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		inline static void sub
		( SmallBinaryFieldBivariatePolynomial & h ,
		  const SmallBinaryFieldBivariatePolynomial & f ,
		  const SmallBinaryFieldBivariatePolynomial & g ) {
			add(h,f,g);
		}

		/**
		 * @brief
		 *            Subtracts an univariate polynomial from the lowest
		 *            Y-coefficient of a bivariate polynomial.
		 *
		 * @details
		 *            This method computes
		 *            \f[
		 *             h(X,Y)=f(X,Y)-s(X).
		 *            \f]
		 *
		 * @param h
		 *            Output polynomial that will contain the difference.
		 *
		 * @param f
		 *            Minuend.
		 *
		 * @param s   Subtrahend.
		 *
		 * @see \link thimble::sub(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldPolynomial&) add()\endlink
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		inline static void sub
		( SmallBinaryFieldBivariatePolynomial & h ,
		  const SmallBinaryFieldBivariatePolynomial & f ,
		  const SmallBinaryFieldPolynomial & s ) {
			add(h,f,s);
		}

		/**
		 * @brief
		 *            Multiplies a bivariate polynomial with a constant
		 *            element of the finite field.
		 *
		 * @details
		 *            For \f$f(X,Y)\in F[X]\f$ and \f$s\in F\f$ this method
		 *            computes
		 *            \f[
		 *             h(X,Y)=s\cdot f(X,Y).
		 *            \f]
		 *
		 * @param h
		 *            The product.
		 *
		 * @param f
		 *            Bivariate factor polynomial
		 *
		 * @param s
		 *            Constant field element.
		 *
		 * @see \link thimble::mul(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,uint32_t) mul()\endlink
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <code>s</code> does not represent a valid element
		 *            of the finite field on which the polynomials
		 *            <code>h</code> and <code>f</code> are based, the
		 *            program runs into undocumented behavior.
		 */
		static void mul
		( SmallBinaryFieldBivariatePolynomial & h ,
		  const SmallBinaryFieldBivariatePolynomial & f ,
		  uint32_t s );

		/**
		 * @brief
		 *            Computes the product of a bivariate polynomial
		 *            with a univariate polynomial.
		 *
		 * @details
		 *            For \f$f(X,Y)\in F[X,Y]\f$ and \f$s(X)\in F[X]\f$ the
		 *            method computes
		 *            \f[
		 *             h(X,Y)=s(X)\cdot g(X,Y).
		 *            \f]
		 *
		 * @param h
		 *            On output, represents the product.
		 *
		 * @param f
		 *            First factor.
		 *
		 * @param s
		 *            Second factor.
		 *
		 * @see \link thimble::mul(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldPolynomial&) mul()\endlink
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		static void mul
		( SmallBinaryFieldBivariatePolynomial & h ,
		  const SmallBinaryFieldBivariatePolynomial & f ,
		  const SmallBinaryFieldPolynomial & s );

		/**
		 * @brief
		 *            Computes the product of two bivariate polynomials.
		 *
		 * @details
		 *            This method computes
		 *            \f[
		 *             h(X,Y)=f(X,Y)\cdot g(X,Y).
		 *            \f]
		 *
		 * @param h
		 *            The product.
		 *
		 * @param f
		 *            First factor.
		 *
		 * @param g
		 *            Second factor.
		 *
		 * @see \link thimble::mul(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&) mul()\endlink
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		static void mul
		( SmallBinaryFieldBivariatePolynomial & h ,
		  const SmallBinaryFieldBivariatePolynomial & f ,
		  const SmallBinaryFieldBivariatePolynomial & g );

		/**
		 * @brief
		 *            Plugs a bivariate polynomial for the Y-value
		 *            of another polynomial and computes the result.
		 *
		 * @details
		 *            This method computes
		 *            \f[
		 *             h(X,Y)=f(X,g(X,Y)).
		 *            \f]
		 *
		 * @param h
		 *            \f$f(X,g(X,Y))\f$.
		 *
		 * @param f
		 *            Bivariate polynomial.
		 *
		 * @param g
		 *            Bivariate polynomial.
		 *
		 * @see \link thimble::evalY(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&) evalY()\endlink
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		static void evalY
		( SmallBinaryFieldBivariatePolynomial & h ,
		  const SmallBinaryFieldBivariatePolynomial & f ,
		  const SmallBinaryFieldBivariatePolynomial & g );

	};

	/**
	 * @brief
	 *            Computes the representation of a bivariate polynomial
	 *            after exchanging its two indeterminate symbols.
	 *
	 * @param g
	 *            The input polynomial \f$g(X,Y)\f$.
	 *
	 * @param f
	 *            The output polynomial which will be
	 *            equals \f$f(X,Y)=g(Y,X)\f$.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	inline void swapXY
	( SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ) {
		SmallBinaryFieldBivariatePolynomial::swapXY(f,g);
	}

	/**
	 * @brief
	 *            Determines the leading term of the highest
	 *            <i>(a,b)</i>-weighted coffient.
	 *
	 * @details
	 *            Write for this polynomial
	 *            \f[
	 *             \sum_{i,j}c_{i,j}X^iY^j.
	 *            \f]
	 *            This function determines <i>(i,j)</i> for which the
	 *            expression <i>a*i+b*j</i> is maximal and then sets
	 *            <code>lt</code> the a representation of
	 *            \f[
	 *             c_{i,j}X^iY^j.
	 *            \f]
	 *            The pair <i>(i,j)</i> will be furthermore returned.
	 *
	 * @param lt
	 *            Reference to which the leading term, according to the
	 *            <i>(a,b)</i>-weighted degree, will be stored.
	 *
	 * @param f
	 *            Input bivariate polynomial from which the leading
	 *            term is computed.
	 *
	 * @param a
	 *            First weight.
	 *
	 * @param b
	 *            Second weight.
	 *
	 * @see deg()
	 *
	 * @return
	 *            The indices <i>(i,j)</i> at which this polynomial's
	 *            cofficient is non-zero and for which the expression
	 *            <i>a*i+b*j</i> is maximal.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	inline std::pair<int,int> leadTerm
	( SmallBinaryFieldBivariatePolynomial & lt ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  int a , int b ) {
		return SmallBinaryFieldBivariatePolynomial::leadTerm(lt,f,a,b);
	}

	/**
	 * @brief
	 *            Essentially, multiplies a bivariate polynomial with
	 *            the factor \f$Y^n\f$.
	 *
	 * @details
	 *            Write
	 *            \f[
	 *             g(X,Y) = \sum_{i,j\geq 0}c_{i,j}X^iY^j.
	 *            \f]
	 *            This method computes
	 *            \f[
	 *             f(X,Y) = Y^n\cdot g(X,Y)
	 *            \f]
	 *            if \f$n\f$ is positive and
	 *            \f[
	 *             f(X,Y) = \sum_{i\geq 0,j\geq -n}c_{i,j}X^iY^{j+n}
	 *            \f]
	 *            if \f$n\f$ is negative.
	 *
	 * @param f
	 *            The output polynomial being the product.
	 *
	 * @param g
	 *            The input factor polynomial.
	 *
	 * @param n
	 *            The exponent of the factor monomial.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	inline void leftShiftY
	( SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ,
	  int n ) {

		SmallBinaryFieldBivariatePolynomial::leftShiftY(f,g,n);
	}

	/**
	 * @brief
	 *            Essentially, computes the quotient by a division
	 *            of a bivariate polynomial divided by \f$Y^n\f$.
	 *
	 * @details
	 *            Write
	 *            \f[
	 *             g(X,Y) = \sum_{i,j\geq 0}c_{i,j}X^iY^j.
	 *            \f]
	 *            This method computes the polynomial
	 *            \f[
	 *             f(X,Y) = \sum_{i\geq 0,j\geq n}c_{i,j}X^iY^{j-n}.
	 *            \f]
	 *            for positive \f$n\f$. If \f$n\f$ is negative, then
	 *            \f[
	 *             f(X,Y) = Y^{-n}\cdot g(X,Y)
	 *            \f]
	 *            is computed.
	 *
	 * @param f
	 *            Output bivariate polynomial the will contain the
	 *            quotient.
	 *
	 * @param g
	 *            Input polynomial divided by \f$Y^n\f$.
	 *
	 * @param n
	 *            The exponent of the divisor monomial.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	inline void rightShiftY
	( SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ,
	  int n ) {

		SmallBinaryFieldBivariatePolynomial::rightShiftY(f,g,n);
	}

	/**
	 * @brief
	 *            Essentially, multiplies a bivariate polynomial with
	 *            the factor \f$X^n\f$.
	 *
	 * @details
	 *            Write
	 *            \f[
	 *             g(X,Y) = \sum_{i,j\geq 0}c_{i,j}X^iY^j.
	 *            \f]
	 *            This method computes
	 *            \f[
	 *             f(X,Y) = X^n\cdot g(X,Y)
	 *            \f]
	 *            if \f$n\f$ is positive and
	 *            \f[
	 *             f(X,Y) = \sum_{i\geq -n,j\geq 0}c_{i,j}X^{i+n}Y^{j}
	 *            \f]
	 *            if \f$n\f$ is negative.
	 *
	 * @param f
	 *            The output polynomial being the product.
	 *
	 * @param g
	 *            The input factor polynomial.
	 *
	 * @param n
	 *            The exponent of the factor monomial.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	inline void leftShiftX
	( SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ,
	  int n ) {

		SmallBinaryFieldBivariatePolynomial::leftShiftX(f,g,n);
	}

	/**
	 * @brief
	 *            Essentially, computes the quotient by a division
	 *            of a bivariate polynomial divided by \f$X^n\f$.
	 *
	 * @details
	 *            Write
	 *            \f[
	 *             g(X,Y) = \sum_{i,j\geq 0}c_{i,j}X^iY^j.
	 *            \f]
	 *            This method computes the polynomial
	 *            \f[
	 *             f(X,Y) = \sum_{i\geq n,j\geq 0}c_{i,j}X^{i-n}Y^{j}.
	 *            \f]
	 *            for positive \f$n\f$. If \f$n\f$ is negative, then
	 *            \f[
	 *             f(X,Y) = X^{-n}\cdot g(X,Y)
	 *            \f]
	 *            is computed.
	 *
	 * @param f
	 *            Output bivariate polynomial the will contain the
	 *            quotient.
	 *
	 * @param g
	 *            Input polynomial divided by \f$X^n\f$.
	 *
	 * @param n
	 *            The exponent of the divisor monomial.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	inline void rightShiftX
	( SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ,
	  int n ) {

		SmallBinaryFieldBivariatePolynomial::rightShiftX(f,g,n);
	}

	/**
	 * @brief
	 *            Computes the sum of two bivariate polynomials.
	 *
	 * @details
	 *            This method computes
	 *            \f[
	 *             h(X,Y)=f(X,Y)+g(X,Y).
	 *            \f]
	 *
	 * @param h
	 *            Output polynomial that will contain the sum.
	 *
	 * @param f
	 *            First summand.
	 *
	 * @param g
	 *            Second summand.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	inline void add
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ) {
		SmallBinaryFieldBivariatePolynomial::add(h,f,g);
	}

	/**
	 * @brief
	 *            Adds an univariate polynomial to the lowest
	 *            Y-coefficient of a bivariate polynomial.
	 *
	 * @details
	 *            This method computes
	 *            \f[
	 *             h(X,Y)=f(X,Y)+s(X).
	 *            \f]
	 *
	 * @param h
	 *            Output polynomial that will contain the sum.
	 *
	 * @param f
	 *            First bivariate summand.
	 *
	 * @param s   Second univariate summand.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	inline void add
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldPolynomial & s ) {
		SmallBinaryFieldBivariatePolynomial::add(h,f,s);
	}

	/**
	 * @brief
	 *            Computes difference between two bivariate polynomials
	 *            which, actually, is equivalent to the computation of
	 *            their sum in binary fields.
	 *
	 * @details
	 *            This method computes
	 *            \f[
	 *             h(X,Y)=f(X,Y)-g(X,Y)
	 *            \f]
	 *            which is equivalent to the computation of their sum.
	 *            Nonetheless, we provide subtraction and negation
	 *            methods explicitly in order to allow the user of
	 *            this library to produce possibly better readable
	 *            code.
	 *
	 * @param h
	 *            Output polynomial that will contain the sum.
	 *
	 * @param f
	 *            Minuend.
	 *
	 * @param g
	 *            Subtrahend.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	inline void sub
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ) {
		SmallBinaryFieldBivariatePolynomial::sub(h,f,g);
	}

	/**
	 * @brief
	 *            Subtracts an univariate polynomial from the lowest
	 *            Y-coefficient of a bivariate polynomial.
	 *
	 * @details
	 *            This method computes
	 *            \f[
	 *             h(X,Y)=f(X,Y)-s(X).
	 *            \f]
	 *
	 * @param h
	 *            Output polynomial that will contain the difference.
	 *
	 * @param f
	 *            Minuend.
	 *
	 * @param s   Subtrahend.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	inline void sub
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldPolynomial & s ) {
		SmallBinaryFieldBivariatePolynomial::sub(h,f,s);
	}

	/**
	 * @brief
	 *            Multiplies a bivariate polynomial with a constant
	 *            element of the finite field.
	 *
	 * @details
	 *            For \f$f(X,Y)\in F[X]\f$ and \f$s\in F\f$ this method
	 *            computes
	 *            \f[
	 *             h(X,Y)=s\cdot f(X,Y).
	 *            \f]
	 *
	 * @param h
	 *            The product.
	 *
	 * @param f
	 *            Bivariate factor polynomial
	 *
	 * @param s
	 *            Constant field element.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If <code>s</code> does not represent a valid element
	 *            of the finite field on which the polynomials
	 *            <code>h</code> and <code>f</code> are based, the
	 *            program runs into undocumented behavior.
	 */
	inline void mul
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  uint32_t s ) {
		SmallBinaryFieldBivariatePolynomial::mul(h,f,s);
	}

	/**
	 * @brief
	 *            Computes the product of a bivariate polynomial
	 *            with a univariate polynomial.
	 *
	 * @details
	 *            For \f$f(X,Y)\in F[X,Y]\f$ and \f$s(X)\in F[X]\f$ the
	 *            method computes
	 *            \f[
	 *             h(X,Y)=s(X)\cdot g(X,Y).
	 *            \f]
	 *
	 * @param h
	 *            On output, represents the product.
	 *
	 * @param f
	 *            First factor.
	 *
	 * @param s
	 *            Second factor.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	inline void mul
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldPolynomial & s ) {
		SmallBinaryFieldBivariatePolynomial::mul(h,f,s);
	}

	/**
	 * @brief
	 *            Computes the product of two bivariate polynomials.
	 *
	 * @details
	 *            This method computes
	 *            \f[
	 *             h(X,Y)=f(X,Y)\cdot g(X,Y).
	 *            \f]
	 *
	 * @param h
	 *            The product.
	 *
	 * @param f
	 *            First factor.
	 *
	 * @param g
	 *            Second factor.
	 *
	 * @see \link thimble::mul(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&) mul()\endlink
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	inline void mul
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ) {
		SmallBinaryFieldBivariatePolynomial::mul(h,f,g);
	}

	/**
	 * @brief
	 *            Plugs a bivariate polynomial for the Y-value
	 *            of another polynomial and computes the result.
	 *
	 * @details
	 *            This method computes
	 *            \f[
	 *             h(X,Y)=f(X,g(X,Y)).
	 *            \f]
	 *
	 * @param h
	 *            \f$f(X,g(X,Y))\f$.
	 *
	 * @param f
	 *            Bivariate polynomial.
	 *
	 * @param g
	 *            Bivariate polynomial.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	inline void evalY
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ) {
		SmallBinaryFieldBivariatePolynomial::evalY(h,f,g);
	}

	/**
	 * @brief
	 *            Prints a text representation of a bivariate polynomial
	 *            to the specified output stream.
	 *
	 * @param out
	 *            The output stream to which a text representation of
	 *            <code>f</code> is written.
	 *
	 * @param f
	 *            The bivariate polynomial of which a text representation
	 *            is written to <code>out</code>.
	 *
	 * @return
	 *            A reference to <code>out</code> after the text
	 *            representation of <code>f</code> has been written
	 *            to <code>out</code>.
	 */
	THIMBLE_DLL std::ostream &operator<<
	( std::ostream & out , const SmallBinaryFieldBivariatePolynomial & f );
}


#endif /* THIMBLE_SMALLBINARYFIELDBIVARIATEPOLYNOMIAL_H_ */
