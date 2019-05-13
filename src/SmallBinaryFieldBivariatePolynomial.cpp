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
 * @file SmallBinaryFieldBivariatePolynomial.cpp
 *
 * @brief
 *            Implementation of a class of which objects represent bivariate
 *            polynomials over small binary galois fields. This class is
 *            provided by the 'SmallBinaryFieldPolynomial.h' header file.
 *
 * @author Benjamin Tams
 */

#include <cstdlib>

#include <thimble/math/numbertheory/SmallBinaryField.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>
#include <thimble/math/numbertheory/SmallBinaryFieldBivariatePolynomial.h>

using namespace std;

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Ensures that \link degreeY\endlink correctly encodes
	 *            the Y-degree of this polynomial.
	 */
	void SmallBinaryFieldBivariatePolynomial::normalize() {

		for ( int i = this->degreeY ; i >= 0 ; i-- ) {
			if ( this->coefficientsY[i].isZero() ) {
				this->degreeY -= 1;
			} else {
				break;
			}
		}
	}

	static void TradMul
	( SmallBinaryFieldPolynomial *c ,
	  const SmallBinaryFieldPolynomial *a  , int m ,
	  const SmallBinaryFieldPolynomial *b , int n ,
	  const SmallBinaryField & gf ) {

		int l = m+n-1;

		for ( int i = 0 ; i < l ; i++ ) {
			c[i].setZero();
		}

		SmallBinaryFieldPolynomial d(gf);

		for ( int i = 0 ; i < m ; i++  ) {

			for ( long j = 0 ; j < n ; j++ ) {
				mul(d,a[i],b[j]);
				add(c[j+i],c[j+i],d);
			}
		}
	}

	/**
	 * @brief
	 *            Computes the product of two bivariate polynomials.
	 *
	 * @details
	 *            see 'SmallBinaryFieldBivariatePolynomial.h'
	 *
	 */
	void SmallBinaryFieldBivariatePolynomial::mulUncheck
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ) {

		SmallBinaryFieldPolynomial *a , *b , *c;
		int m , n;

		if ( f.degY() < g.degY() ) {
			m = f.degY();
			n = g.degY();
			a = f.coefficientsY;
			b = g.coefficientsY;
		} else {
			m = g.degY();
			n = f.degY();
			a = g.coefficientsY;
			b = f.coefficientsY;
		}

		h.ensureCapacity(m+n+1);
		h.degreeY = m+n;

		c = h.coefficientsY;

		TradMul(c,a,m+1,b,n+1,h.gfPtr[0]);

		h.normalize();
	}

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
	void SmallBinaryFieldBivariatePolynomial::ensureCapacity( int newCapacity ) {

		// Check whether capacity must be increased
		if ( newCapacity > this->capacity ) {

            if ( this->coefficientsY == NULL ) {
                this->coefficientsY = (SmallBinaryFieldPolynomial*)malloc
                        ( newCapacity * sizeof(SmallBinaryFieldPolynomial) );
            } else {
                this->coefficientsY = (SmallBinaryFieldPolynomial*)realloc
                        ( this->coefficientsY , newCapacity * sizeof(SmallBinaryFieldPolynomial));
            }

            if ( this->coefficientsY == NULL ) {
                cerr << "SmallBinaryFieldBivariatePolynomial::ensureCapacity "
                     << "out of memory." << endl;
                exit(EXIT_FAILURE);
            }

            for ( int i = this->capacity ; i < newCapacity ; i++ ) {
                new (this->coefficientsY+i) SmallBinaryFieldPolynomial(this->gfPtr[0]);
            }

            this->capacity = newCapacity;
		}
	}

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
	void SmallBinaryFieldBivariatePolynomial::assign
	( const SmallBinaryFieldBivariatePolynomial & f ) {

		if ( this == &f ) {
			return;
		}

		if ( f.isZero() ) {
			setZero();
			return;
		}

		ensureCapacity(f.degY()+1);

		this->gfPtr = f.gfPtr;

		for ( int i = 0 ; i <= f.degY() ; i++ ) {
			this->coefficientsY[i].assign(f.coefficientsY[i]);
		}

		this->degreeY = f.degreeY;
	}

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
    SmallBinaryFieldBivariatePolynomial::SmallBinaryFieldBivariatePolynomial
    ( const SmallBinaryFieldPolynomial & c0 ) {

        this->gfPtr = &c0.getField();

        this->coefficientsY = (SmallBinaryFieldPolynomial*)malloc
                ( sizeof(SmallBinaryFieldPolynomial) );
        if ( this->coefficientsY == NULL ) {
            cerr << "SmallBinaryFieldBivariatePolynomial: "
                 << "out of memory." << endl;
            exit(EXIT_FAILURE);
        }

        new (this->coefficientsY) SmallBinaryFieldPolynomial(c0);

        this->degreeY = 0;
        this->capacity = 1;
        normalize();
    }

    /**
     * @brief
     *            Destructor.
     *
     * @details
     *            Frees all the memory that has been allocated to
     *            represent this bivariate polynomial.
     */
    SmallBinaryFieldBivariatePolynomial::~SmallBinaryFieldBivariatePolynomial() {

        for ( int j = 0 ; j < this->capacity ; j++ ) {
            (this->coefficientsY+j)->~SmallBinaryFieldPolynomial();
        }

        free(this->coefficientsY);
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
	std::pair<int,int> SmallBinaryFieldBivariatePolynomial::deg
	( int a , int b ) const {

		bool chosen = false;
		std::pair<long,long> d(-1,-1);

		int dy = degY();

		for ( int j = 0 ; j <= dy ; j++ ) {

			if ( this->coefficientsY[j].isZero() ) {
				continue;
			}

			// Now, 'coefficients[j]' is non-zero.

			int i;

			// If 'a' is negative, ...
			if ( a < 0 ) {
				// ... search the minimal 'i' because this will maximize
				// 'a*i+b*j'.
				for ( i = 0 ; this->coefficientsY[j].getCoeff(i) == 0 ; i++ ) {
				}
			} else {
				// Otherwise, if 'a' is positive, the maximal 'i', which is
				// the X-degree, will maximize the expression 'a*i+b*j'.
				i = this->coefficientsY[j].deg();
				if ( i < 0 ) {
					i = 0;
				}
			}

			// Check if the new (a,b)-degree candidate is larger and if it is
			// or if it has not yet been set, ...
			if ( !chosen || a*d.first+b*d.second < a*i+b*j) {
				// ... update.
				d.first = i;
				d.second = j;
				chosen = true;
			}
		}

		return d;
	}

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
	uint32_t SmallBinaryFieldBivariatePolynomial::getCoeff
	( int i , int j ) const {

		if ( i < 0 || j < 0 ) {
			return 0;
		}

		if ( j > this->degreeY ) {
			return 0;
		}

		return this->coefficientsY[j].getCoeff(i);
	}

	/**
	 * @brief
	 *           Determine whether this bivariate polynomial is equals to the
	 *           given bivariate polynomial.
	 *
	 * @details
	 *           see 'SmallBinaryFieldBivariatePolynomial.h'
	 */
	bool SmallBinaryFieldBivariatePolynomial::equals
	( const SmallBinaryFieldBivariatePolynomial & f ) const {

		// If 'f' is of the same reference as this polynomial
		// the polynomials must be equal
		if ( this == &f ) {
			return true;
		}

		// Check whether the two polynomials are defined over the same
		// finite field.
		if ( this->gfPtr != f.gfPtr ) {
			cerr << "SmallBinaryFieldBivariatePolynomial::equals: When two "
				 << "polynomials are compared they must have coefficients in "
				 << "the same field."
				 << endl;
			exit(EXIT_FAILURE);
		}

		// If the degrees differ the polynomial can not be equal
		int n = this->degreeY;
		if ( n != f.degreeY ) {
			return false;
		}

		// Compare the relevant Y-coefficients.
		for ( int i = 0 ; i <= n ; i++ ) {
			if ( !(this->coefficientsY[i].equals(f.coefficientsY[i])) ) {
				return false;
			}
		}

		return true;
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
	void SmallBinaryFieldBivariatePolynomial::setCoeffY
	( int i , const SmallBinaryFieldPolynomial & c ) {

		// Check whether a valid index is accessed ...
		if ( i < 0 ) {
			// ... if not, verbose and exit.
			cerr << "SmallBinaryFieldBivariatePolynomial::setCoeff: Accessed "
				 << "coefficient must have positive index." << endl;
			exit(EXIT_FAILURE);
		}

		if ( this->gfPtr->getDefiningPolynomial().rep !=
			 c.gfPtr->getDefiningPolynomial().rep ) {

			cerr << "SmallBinaryFieldBivariatePolynomial::setCoeff: "
				 << "conflicting finite field context." << endl;
			exit(EXIT_FAILURE);
		}

		// Save the Y-degree of the polynomial
		int d = this->degreeY;

		if ( d < i ) {

			// If a coefficient beyond relevant indices
			// is accessed, ensure sufficient capacity, ...
			ensureCapacity(i+1);

			// ..., ensure the offset to be zero, and ...
			for ( int j = d+1 ; j < i ; j++ ) {
				this->coefficientsY[j].setZero();
			}

			// ..., and update the Y-degree
			this->degreeY = i;
		}

		// Update the accessed coefficient.
		this->coefficientsY[i].assign(c);

		// The degree might be wrong if accessed index is
		// higher than the old degree of the polynomial
		if ( d <= i ) {
			normalize();
		}
	}

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
	void SmallBinaryFieldBivariatePolynomial::setCoeff
	( int i , int j , uint32_t c ) {

		if ( i < 0 || j < 0 ) {
			cerr << "SmallBinaryFieldBivariatePolynomial::setCoeff: "
			     << "Accessed coefficient's indices must be positive." << endl;
			exit(EXIT_FAILURE);
		}

		if ( c >= this->gfPtr[0].getCardinality() ) {
			cerr << "SmallBinaryFieldBivariatePolynomial::setCoeff: "
				 << "specified coefficient does not represent a valid "
				 << "element in the finite field context." << endl;
			exit(EXIT_FAILURE);
		}

		int d = this->degreeY;

		if ( d < j ) {

			ensureCapacity(j+1);

			for ( int l = d+1 ; l <= j ; l++ ) {
				this->coefficientsY[l].setZero();
			}

			this->degreeY = j;
		}


		this->coefficientsY[j].setCoeff(i,c);

		this->normalize();
	}

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
	uint32_t SmallBinaryFieldBivariatePolynomial::eval
	( uint32_t x , uint32_t y ) const {

	    int dy = degY();

	    if ( dy >= 0 ) {

	    	SmallBinaryFieldPolynomial f(this->coefficientsY[dy]);

	    	// Horner rule like evaluation at 'y' ...
	    	for ( int i = dy-1 ; i >= 0 ; i-- ) {
	    		SmallBinaryFieldPolynomial::mul(f,f,y);
	    		SmallBinaryFieldPolynomial::add(f,f,this->coefficientsY[i]);
	    	}

	    	// ... and finally at 'x'.
	    	return f.eval(x);
	    }

	    return 0;
	}

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
	void SmallBinaryFieldBivariatePolynomial::swapXY
	( SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ) {

		if ( &f == &g ) {
			SmallBinaryFieldBivariatePolynomial tg(g);
			swapXY(f,tg);
			return;
		}

		f.setZero();

		uint32_t c;

		int dy = g.degY();

		for ( int j = dy ; j >= 0 ; j-- ) {
			int dx = g.degX(j);
			for ( int i = dx ; i >= 0 ; i-- ) {
				c = g.getCoeff(i,j);
				f.setCoeff(j,i,c);
			}
		}
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
	pair<int,int> SmallBinaryFieldBivariatePolynomial::leadTerm
	( SmallBinaryFieldBivariatePolynomial & lt ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  int a , int b ) {

		  pair<int,int> d = f.deg(a,b);

		  if ( d.second >= 0 && d.first >= 0) {
			  uint32_t c = f.getCoeff(d.first,d.second);
			  lt.setZero();
			  lt.setCoeff(d.first,d.second,c);
		  } else {
			  lt.setZero();
		  }

		  return d;
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
	 * @see \link thimble::leftShiftY(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,int) leftShiftY()\endlink
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void SmallBinaryFieldBivariatePolynomial::leftShiftY
	( SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ,
	  int n ) {

		if ( n < 0 ) {
			rightShiftY(f,g,-n);
			return;
		}

		if ( !g.isZero() && n > 0 ) {

			int d = g.degreeY;

			f.ensureCapacity(d+n+1);
			f.degreeY = d+n;

			SmallBinaryFieldPolynomial *a , *b;
			a = f.coefficientsY;
			b = g.coefficientsY;

			for ( int i = d ; i >= 0 ; i-- ) {
				a[i+n].assign(b[i]);
			}

			a = f.coefficientsY;
			for ( int i = 0 ; i < n ; i++ ) {
				a[i].setZero();
			}
		} else {
			f.assign(g);
		}
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
	 * @see \link thimble::rightShiftY(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,int) rightShiftY()\endlink
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void SmallBinaryFieldBivariatePolynomial::rightShiftY
	( SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ,
	  int n ) {

		if ( n < 0 ) {
			leftShiftY(f,g,n);
			return;
		}

		f.assign(g);

		int d = f.degY();

		SmallBinaryFieldPolynomial *a = f.coefficientsY;

		for ( int i = 0 ; i+n <= d  ; i++ ) {
			a[i] = a[i+n];
		}

		if ( (f.degreeY -= n) < 0 ) {
			f.degreeY = -1;
		}
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
	 * @see \link thimble::leftShiftX(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,int) leftShiftX()\endlink
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void SmallBinaryFieldBivariatePolynomial::leftShiftX
	( SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ,
	  int n ) {

		f.assign(g);

		int dy = f.degY();

		SmallBinaryFieldPolynomial *a = f.coefficientsY;
		for ( int j = 0 ; j <= dy ; j++ , a++ ) {
			leftShift(*a,*a,n);
		}

		f.normalize();
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
	 * @see \link thimble::rightShiftX(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,int) rightShiftX()\endlink
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void SmallBinaryFieldBivariatePolynomial::rightShiftX
	( SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ,
	  int n ) {

		f.assign(g);

		int dy = f.degY();

		SmallBinaryFieldPolynomial *a = f.coefficientsY;
		for ( int j = 0 ; j <= dy ; j++ , a++ ) {
			rightShift(*a,*a,n);
		}

		f.normalize();
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
	void SmallBinaryFieldBivariatePolynomial::add
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ) {

		SmallBinaryFieldPolynomial *a , *b , *c;
		int m , n;

		// Obtain the higher Y-degree 'n' and the lower
		// Y-degree of the summands.
		if ( f.degY() < g.degY() ) {
			m = f.degY();
			n = g.degY();
		} else {
			m = g.degY();
			n = f.degY();
		}

		// Ensure the output polynomial can store the sum.
		h.ensureCapacity(n+1);

		// Ensure that 'a' contains the coefficients of the summand
		// that is of the higher degree 'n'; the coefficients of the
		// other polynomial of the lower degree 'm' are accessed by 'b'.
		if ( f.degY() < g.degY() ) {
			b = f.coefficientsY;
			c = g.coefficientsY;
		} else {
			b = g.coefficientsY;
			c = f.coefficientsY;
		}

		// Temporarily, set the maximal possible degree.
		h.degreeY = n;

		// 'a' accesses to coefficient of the output polynomial
		a = h.coefficientsY;

		// The first 'm+1' coefficients of the output polynomial
		// corresponds to the sum of the first respective 'm+1' elements
		// of the two input polynomials. Each of these sums, in turn,
		// correspond to a binary exclusive or operation.
		for ( int i = 0 ; i <= m ; i++ ) {
			thimble::add(a[i],b[i],c[i]);
		}

		// The remaining coefficients simply are a copy of the corresponding
		// coefficient of the higher degree polynomial.
		if ( a != c && n > m ) {
			for ( int i = m+1 ; i <= n ; i++ ) {
				a[i].assign(c[i]);
			}
		}

		// The degree might be wrong, e.g., if the highest index coefficient
		// summed to zero. Therefore, ensure the Y-degree is correct.
		h.normalize();
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
	 * @see \link thimble::add(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldPolynomial&) add()\endlink
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void SmallBinaryFieldBivariatePolynomial::add
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldPolynomial & g ) {

		h.assign(f);

		if ( h.isZero() ) {
			h.ensureCapacity(1);
			h.coefficientsY[0].assign(g);
			h.degreeY = 0;
		} else {
			SmallBinaryFieldPolynomial::add
				(h.coefficientsY[0],h.coefficientsY[0],g);
		}

		h.normalize();
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
	 * @see \link thimble::mul(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial,const SmallBinaryFieldPolynomial) mul()\endlink
	 *
	 * @warning
	 *            If <code>s</code> does not represent a valid element
	 *            of the finite field on which the polynomials
	 *            <code>h</code> and <code>f</code> are based, the
	 *            program runs into undocumented behavior.
	 */
	void SmallBinaryFieldBivariatePolynomial::mul
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  uint32_t s ) {

		h.assign(f);

		for ( int i = 0 ; i <= h.degreeY ; i++ ) {
			SmallBinaryFieldPolynomial::mul
				(h.coefficientsY[i],h.coefficientsY[i],s);
		}

		h.normalize();
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
	 * @see \link thimble::mul(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldPolynomial&) mul()\endlink
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void SmallBinaryFieldBivariatePolynomial::mul
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldPolynomial & s ) {

		int dy = f.degY();

		h.ensureCapacity(dy+1);
		h.degreeY = dy;

		for ( int j = 0 ; j <= dy ; j++ ) {
			SmallBinaryFieldPolynomial::mul
				(h.coefficientsY[j],f.coefficientsY[j],s);
		}

		h.normalize();
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
	void SmallBinaryFieldBivariatePolynomial::mul
	( SmallBinaryFieldBivariatePolynomial & h ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & g ) {

		if ( h.gfPtr != f.gfPtr || h.gfPtr != g.gfPtr ) {
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
			SmallBinaryFieldBivariatePolynomial tf(f) , tg(g);
			mulUncheck(h,tf,tg);
		} else if ( &h == &f ) {
			SmallBinaryFieldBivariatePolynomial tf(f);
			mulUncheck(h,tf,g);
		} else if ( &h == &g ) {
			SmallBinaryFieldBivariatePolynomial tg(g);
			mulUncheck(h,f,tg);
		} else {
			mulUncheck(h,f,g);
		}
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
	 * @see \link thimble::evalY(SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&,const SmallBinaryFieldBivariatePolynomial&) evalY()\endlink
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void SmallBinaryFieldBivariatePolynomial::evalY
	( SmallBinaryFieldBivariatePolynomial & g ,
	  const SmallBinaryFieldBivariatePolynomial & f ,
	  const SmallBinaryFieldBivariatePolynomial & y ) {

		if ( &g == &f ) {
			SmallBinaryFieldBivariatePolynomial tf(f);
			evalY(g,tf,y);
			return;
	    } else if ( &g == &y ) {
	    	SmallBinaryFieldBivariatePolynomial ty(y);
	    	evalY(g,f,ty);
	    	return;
	    }

	    int dy = f.degY();

	    g.setZero();

	    if ( dy >= 0 ) {

	    	g.setCoeffY(0,f.coefficientsY[dy]);

	    	for ( int j = dy-1 ; j >= 0 ; j-- ) {
	    		mul(g,g,y);
	    		add(g,g,f.coefficientsY[j]);
	    	}
	    }
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
	ostream &operator<<
	( ostream & out , const SmallBinaryFieldBivariatePolynomial & f ) {

	    int d = f.degY();

	    if ( d >= 0 ) {

	    	// If the degree is larger than 0, then at least
	    	// one element is non-zero causing a non-empty printout

	    	for ( int i = 0 ; i <= d ; i++ ) {

	    		if ( !f.getCoeffY(i).isZero() ) {

	    			out << "(0x" << std::hex << f.getCoeffY(i) << ")";
	    			if ( i > 0 ) {
	    				out << " * Y";
	    				if ( i > 1 )
	    					out << "^" << std::dec << i;
	    			}

	    			if ( i < d )
	    				out << " + ";
	    		}
	      }

	    } else {

	    	// If the polynomial is zero, print a hexadecimal zero.
	    	out << "0x0";
	    }

	    out << std::dec;

		return out;

	}
}



